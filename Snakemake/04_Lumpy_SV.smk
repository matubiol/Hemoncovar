"""
PIPELINE: Lumpy SV Detection & Clinical Annotation
DESCRIPTION:
This Snakemake workflow performs structural variant detection using Lumpy,
followed by gene annotation of both breakpoints (04.1_Gene_Ann_SV.py) and
clinical annotation using the CIViC database (04.2_CIVic_Ann_SV.R).
There is no Snakemake wrapper for Lumpy: it runs in a pinned conda environment.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Extract split reads and discordant pairs.
rule extract_reads:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai"
    output:
        split = temp(f"{LUMPY_DIR}/01_SV_Calling/{{sample}}_splitters.bam"),
        disc  = temp(f"{LUMPY_DIR}/01_SV_Calling/{{sample}}_discordants.bam")
    log:
        f"{LUMPY_DIR}/01_SV_Calling/Logs/{{sample}}_extract_reads.log"
    conda:
        ENVS['lumpy']
    threads: 8
    resources:
        mem_mb=20000
    shell:
        """
        set -o pipefail
        # Split reads: take ALL reads. Split reads frequently belong to
        # properly-paired templates, so -F 1294 must NOT be applied here.
        samtools view -h {input.bam} 2> {log} \
            | extractSplitReads_BwaMem -i stdin 2>> {log} \
            | samtools view -Sb - \
            | samtools sort -@ {threads} -T {resources.tmpdir}/{wildcards.sample}_split -o {output.split} 2>> {log}

        # Discordant pairs: exclude proper pairs, unmapped reads/mates, secondary and duplicates
        samtools view -b -F 1294 {input.bam} 2>> {log} \
            | samtools sort -@ {threads} -T {resources.tmpdir}/{wildcards.sample}_disc -o {output.disc} 2>> {log}
        """

# Run LUMPY structural variant caller.
rule run_lumpy:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        split = f"{LUMPY_DIR}/01_SV_Calling/{{sample}}_splitters.bam",
        disc = f"{LUMPY_DIR}/01_SV_Calling/{{sample}}_discordants.bam",
        ref = GENOME
    output:
        vcf = f"{LUMPY_DIR}/01_SV_Calling/{{sample}}.vcf"
    log:
        f"{LUMPY_DIR}/01_SV_Calling/Logs/{{sample}}.log",
    benchmark:
        f"{LUMPY_DIR}/01_SV_Calling/Logs/{{sample}}.benchmark"
    conda:
        ENVS['lumpy']
    resources:
        mem_mb=20000
    shell:
        """
        lumpyexpress \
            -B {input.bam} \
            -R {input.ref} \
            -S {input.split} \
            -D {input.disc} \
            -o {output.vcf} 2> {log}
        """

# Gene annotation of both breakpoints
rule ann_genes_lumpy:
    input:
        vcf = f"{LUMPY_DIR}/01_SV_Calling/{{sample}}.vcf",
        gtf = REFS['gtf']
    params:
        min_support = config.get("sv_min_support", 3),
        pass_only = config.get("sv_pass_only", True),
        max_spanned_genes = config.get("sv_max_spanned_genes", 50)
    output:
        tsv = temp(f"{LUMPY_DIR}/02_Annotated/{{sample}}_ann_gene.tsv")
    log:
        f"{LUMPY_DIR}/02_Annotated/Logs/{{sample}}_ann_gene.log"
    conda:
        ENVS['python']
    script:
        SRC['SV_gene_ann']

# Clinical annotation using CIViC database
rule ann_civic_lumpy:
    input:
        unpack(civic_inputs(LUMPY_DIR))
    output:
        xlsx = f"{LUMPY_DIR}/02_Annotated/{{sample}}_ann_civic.xlsx"
    log:
        f"{LUMPY_DIR}/02_Annotated/Logs/{{sample}}_ann_civic.log"
    conda:
        ENVS['r_civic']
    script:
        SRC['SV_CIVIC_ann']
