"""
PIPELINE: Manta SV Detection & Clinical Annotation
DESCRIPTION:
This Snakemake workflow performs tumour-only structural variant detection
using Manta (exome/targeted mode), followed by gene annotation of both
breakpoints (04.1_Gene_Ann_SV.py) and clinical annotation using the CIViC
database (04.2_CIVic_Ann_SV.R).
The Manta wrapper is NOT used: it always passes the BAMs as --normalBam
(germline diploid model) and returns diploidSV.vcf.gz, whereas tumour-only
somatic calling needs --tumorBam and tumorSV.vcf.gz. Manta runs instead in a
pinned conda environment (manta 1.6.0, same version as the wrapper).
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Configure and run Manta
rule run_manta:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        db = GENOME,
        fai = GENOME_FAI
    params:
        dir_name = f"{MANTA_DIR}/01_SV_Calling/{{sample}}",
        script_dir = f"{MANTA_DIR}/01_SV_Calling/{{sample}}/runWorkflow.py"
    output:
        vcf = f"{MANTA_DIR}/01_SV_Calling/{{sample}}/results/variants/tumorSV.vcf.gz"
    log:
        f"{MANTA_DIR}/01_SV_Calling/Logs/{{sample}}.log",
    benchmark:
        f"{MANTA_DIR}/01_SV_Calling/Logs/{{sample}}.benchmark"
    conda:
        ENVS['manta']
    threads: 16
    resources:
        mem_mb=32000
    shell:
        """
        # Remove leftovers of a previous (failed) run: configManta.py refuses
        # to overwrite an existing runWorkflow.py
        rm -rf {params.script_dir} {params.script_dir}.config.pickle {params.dir_name}/workspace

        # Configure Manta for tumour-only SV calling on targeted data
        configManta.py \
            --tumorBam {input.bam} \
            --referenceFasta {input.db} \
            --exome \
            --runDir {params.dir_name} 2> {log}

        # Run Manta workflow
        python2 {params.script_dir} \
            -m local \
            -j {threads} 2>> {log}
        """

# Gene annotation of both breakpoints
rule ann_genes_manta:
    input:
        vcf = f"{MANTA_DIR}/01_SV_Calling/{{sample}}/results/variants/tumorSV.vcf.gz",
        gtf = REFS['gtf']
    params:
        min_support = config.get("sv_min_support", 3),
        pass_only = config.get("sv_pass_only", True),
        max_spanned_genes = config.get("sv_max_spanned_genes", 50)
    output:
        tsv = temp(f"{MANTA_DIR}/02_Annotated/{{sample}}_ann_gene.tsv")
    log:
        f"{MANTA_DIR}/02_Annotated/Logs/{{sample}}_ann_gene.log"
    conda:
        ENVS['python']
    script:
        SRC['SV_gene_ann']

# Clinical annotation using CIViC database
rule ann_civic_manta:
    input:
        unpack(civic_inputs(MANTA_DIR))
    output:
        xlsx = f"{MANTA_DIR}/02_Annotated/{{sample}}_ann_civic.xlsx"
    log:
        f"{MANTA_DIR}/02_Annotated/Logs/{{sample}}_ann_civic.log"
    conda:
        ENVS['r_civic']
    script:
        SRC['SV_CIVIC_ann']
