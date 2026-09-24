"""
PIPELINE: GATK Copy Number Variation Analysis
DESCRIPTION:
Somatic CNV workflow (GATK4, tumour-only against a panel of normals).

  Panel of normals (build mode, written to PON_DIR/gatk):
    targets.preprocessed.interval_list (padded, no binning; chrX/chrY excluded
    by default because the PoN usually mixes sexes), panel_of_normals.hdf5
    built from ALL normal samples.
  Per tumour sample:
    read counts -> denoising -> allelic counts at common SNPs within the
    targets -> joint copy-ratio / allele-fraction segmentation (ModelSegments,
    which also captures copy-neutral LOH) -> segment calling -> plots ->
    gene/cytoband annotation.
  PoN QC (build mode, pon_qc): every normal is denoised, segmented and called
    against a PoN built without it.

Wrappers (GATK 4.6.2.0): CollectReadCounts, CollectAllelicCounts,
DenoiseReadCounts, ModelSegments, CallCopyRatioSegments. Pinned conda
environments (same GATK version): PreprocessIntervals,
CreateReadCountPanelOfNormals and the plotting tools (need R).
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
GATK_PON = f"{PON_DIR}/gatk"
CNV_INTERVALS = f"{GATK_PON}/targets.preprocessed.interval_list"
NEUTRAL_LOW, NEUTRAL_HIGH = CNV_CFG.get("neutral_bounds", [0.9, 1.1])
CALL_EXTRA = (f"--neutral-segment-copy-ratio-lower-bound {NEUTRAL_LOW} "
              f"--neutral-segment-copy-ratio-upper-bound {NEUTRAL_HIGH}")

# =============================================================================
# PANEL OF NORMALS (build mode)
# =============================================================================
if PON_MODE == "build":
    # Padded target intervals (no binning for capture data) - no wrapper
    rule process_int:
        input:
            bed = INTERVALS,
            ref = GENOME,
            dict = GENOME_DICT
        params:
            padding = CNV_CFG.get("interval_padding", 250),
            exclude = " ".join(f"-XL {c}" for c in SEX_CONTIGS) if CNV_CFG.get("gatk_exclude_sex_chroms", True) else ""
        output:
            CNV_INTERVALS
        log:
            f"{GATK_PON}/Logs/preprocess_intervals.log"
        conda:
            ENVS['gatk']
        shell:
            """
            gatk PreprocessIntervals \
                -L {input.bed} \
                {params.exclude} \
                -R {input.ref} \
                --padding {params.padding} \
                --bin-length 0 \
                -imr OVERLAPPING_ONLY \
                -O {output} 2> {log}
            """

    # Panel of normals from ALL normal samples - no wrapper
    rule create_pon_gatk:
        input:
            counts = expand(f"{GATK_CNV_DIR}/00_ReadCounts/{{normal}}.counts.hdf5", normal=NORMALS)
        params:
            inputs = lambda wildcards, input: " ".join(f"-I {c}" for c in input.counts)
        output:
            f"{GATK_PON}/panel_of_normals.hdf5"
        log:
            f"{GATK_PON}/Logs/create_pon.log"
        conda:
            ENVS['gatk']
        resources:
            mem_mb = 16000
        shell:
            """
            gatk --java-options "-Xmx12g" CreateReadCountPanelOfNormals \
                {params.inputs} \
                -O {output} 2> {log}
            """

# =============================================================================
# PER-SAMPLE RULES (normals and tumours)
# =============================================================================
rule read_counts:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        intervals = CNV_INTERVALS
    output:
        counts = f"{GATK_CNV_DIR}/00_ReadCounts/{{sample}}.counts.hdf5"
    params:
        mergingRule = "OVERLAPPING_ONLY"
    log:
        f"{GATK_CNV_DIR}/00_ReadCounts/Logs/{{sample}}.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/collectreadcounts"

# =============================================================================
# TUMOUR SAMPLES
# =============================================================================
# Allelic counts at common SNPs within the (padded) targets
rule allelic_counts:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        ref = GENOME,
        dict = GENOME_DICT,
        intervals = REFS['common_snps'],
        targets = CNV_INTERVALS
    output:
        counts = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.allelicCounts.tsv"
    params:
        extra = lambda wildcards, input: f"-L {input.targets} --interval-set-rule INTERSECTION"
    log:
        f"{GATK_CNV_DIR}/02_Segmentation/Logs/{{sample}}_allelic_counts.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/collectalleliccounts"

# Denoise read counts using the panel of normals
rule denoise_counts:
    input:
        hdf5 = f"{GATK_CNV_DIR}/00_ReadCounts/{{sample}}.counts.hdf5",
        pon = f"{GATK_PON}/panel_of_normals.hdf5",
        pon_check = PON_CHECK
    output:
        std_copy_ratio = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.standardizedCR.tsv",
        denoised_copy_ratio = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoisedCR.tsv"
    log:
        f"{GATK_CNV_DIR}/01_Denoising/Logs/{{sample}}_denoise.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/denoisereadcounts"

# Plot denoised copy ratios - no wrapper (needs R)
rule plot_denoised:
    input:
        scr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.standardizedCR.tsv",
        dcr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoisedCR.tsv",
        dic = GENOME_DICT
    output:
        f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoised.png",
        temp([
            f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.{suffix}"
            for suffix in [
                "deltaMAD.txt", "denoisedMAD.txt",
                "scaledDeltaMAD.txt", "standardizedMAD.txt"
            ]
        ])
    log:
        f"{GATK_CNV_DIR}/01_Denoising/Logs/{{sample}}_plot.log"
    conda:
        ENVS['gatk_plots']
    shell:
        """
        gatk PlotDenoisedCopyRatios \
            --standardized-copy-ratios {input.scr} \
            --denoised-copy-ratios {input.dcr} \
            --sequence-dictionary {input.dic} \
            --minimum-contig-length 46709983 \
            --output-prefix {wildcards.sample} \
            --output {GATK_CNV_DIR}/01_Denoising 2> {log}
        """

# Joint copy-ratio / allele-fraction segmentation
rule model_segments:
    input:
        denoised_copy_ratios = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoisedCR.tsv",
        allelic_counts = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.allelicCounts.tsv"
    output:
        cr_seg = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.cr.seg",
        cr_igv = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.cr.igv.seg",
        af_igv = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.af.igv.seg",
        hets = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.hets.tsv",
        model_final = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.modelFinal.seg"
    log:
        f"{GATK_CNV_DIR}/02_Segmentation/Logs/{{sample}}_model_segments.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/modelsegments"

# Plot modeled segments - no wrapper (needs R)
rule plot_segmented:
    input:
        dcr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoisedCR.tsv",
        ac = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.hets.tsv",
        seg = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.modelFinal.seg",
        dic = GENOME_DICT
    output:
        f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.modeled.png"
    log:
        f"{GATK_CNV_DIR}/02_Segmentation/Logs/{{sample}}_plot.log"
    conda:
        ENVS['gatk_plots']
    shell:
        """
        gatk PlotModeledSegments \
            --denoised-copy-ratios {input.dcr} \
            --allelic-counts {input.ac} \
            --segments {input.seg} \
            --sequence-dictionary {input.dic} \
            --minimum-contig-length 46709983 \
            --output-prefix {wildcards.sample} \
            --output {GATK_CNV_DIR}/02_Segmentation 2> {log}
        """

# Call copy-ratio segments (gain / loss / neutral)
rule call_cnvs:
    input:
        copy_ratio_seg = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.cr.seg"
    output:
        called_copy_ratio_seg = f"{GATK_CNV_DIR}/03_CNV_Calling/{{sample}}.called.seg"
    params:
        extra = CALL_EXTRA
    log:
        f"{GATK_CNV_DIR}/03_CNV_Calling/Logs/{{sample}}_call.log"
    resources:
        mem_mb = 4000
    wrapper:
        f"{WRAPPER}/bio/gatk/callcopyratiosegments"

# Gene (GTF) and cytoband annotation of all segments
rule annotate_coord_gatk:
    input:
        call = f"{GATK_CNV_DIR}/03_CNV_Calling/{{sample}}.called.seg",
        cyto = REFS['cyto'],
        gtf = REFS['gtf']
    params:
        format = "gatk",
        min_probes = 0,
        sample = "{sample}",
        sex = lambda wildcards: SAMPLE_INFO[wildcards.sample]["sex"]
    output:
        tsv = f"{GATK_CNV_DIR}/04_Annotated/{{sample}}.ann.tsv"
    log:
        f"{GATK_CNV_DIR}/04_Annotated/Logs/{{sample}}.log"
    conda:
        ENVS['python']
    script:
        SRC['CNV_gene_ann']

# =============================================================================
# PoN QC: leave-one-out evaluation of the normals (build mode)
# =============================================================================
if PON_QC:
    rule gatk_loo_pon:
        input:
            counts = lambda wildcards: [f"{GATK_CNV_DIR}/00_ReadCounts/{n}.counts.hdf5"
                                        for n in NORMALS if n != wildcards.normal]
        params:
            inputs = lambda wildcards, input: " ".join(f"-I {c}" for c in input.counts)
        output:
            temp(f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.loo_pon.hdf5")
        log:
            f"{GATK_CNV_DIR}/PoN_QC/Logs/{{normal}}_pon.log"
        conda:
            ENVS['gatk']
        resources:
            mem_mb = 16000
        shell:
            """
            gatk --java-options "-Xmx12g" CreateReadCountPanelOfNormals \
                {params.inputs} \
                -O {output} 2> {log}
            """

    rule gatk_loo_denoise:
        input:
            hdf5 = f"{GATK_CNV_DIR}/00_ReadCounts/{{normal}}.counts.hdf5",
            pon = f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.loo_pon.hdf5"
        output:
            std_copy_ratio = temp(f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.standardizedCR.tsv"),
            denoised_copy_ratio = f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.denoisedCR.tsv"
        log:
            f"{GATK_CNV_DIR}/PoN_QC/Logs/{{normal}}_denoise.log"
        resources:
            mem_mb = 8000
        wrapper:
            f"{WRAPPER}/bio/gatk/denoisereadcounts"

    rule gatk_loo_segments:
        input:
            denoised_copy_ratios = f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.denoisedCR.tsv"
        output:
            cr_seg = temp(f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.cr.seg")
        log:
            f"{GATK_CNV_DIR}/PoN_QC/Logs/{{normal}}_segments.log"
        resources:
            mem_mb = 8000
        wrapper:
            f"{WRAPPER}/bio/gatk/modelsegments"

    rule gatk_loo_call:
        input:
            copy_ratio_seg = f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.cr.seg"
        output:
            called_copy_ratio_seg = f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.called.seg"
        params:
            extra = CALL_EXTRA
        log:
            f"{GATK_CNV_DIR}/PoN_QC/Logs/{{normal}}_call.log"
        resources:
            mem_mb = 4000
        wrapper:
            f"{WRAPPER}/bio/gatk/callcopyratiosegments"
