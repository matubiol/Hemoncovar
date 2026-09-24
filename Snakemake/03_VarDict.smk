"""
PIPELINE: VarDict Small Variation Analysis
DESCRIPTION:
Somatic variant calling with VarDict (tumour-only, Snakemake wrapper:
vardict-java | teststrandbias.R | var2vcf_valid.pl). The raw VCF keeps every
call with its FILTER column filled; normalization, annotation and reporting
(PASS calls by default) are done in 03_SmallVariants_Common.smk.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
rule vardict:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        reference = GENOME,
        fai = GENOME_FAI,
        regions = INTERVALS
    output:
        vcf = f"{VAR_DIR}/01_Variants/{{sample}}.raw.vcf"
    params:
        sample_name = "{sample}",
        allele_frequency_threshold = "0.01",
        bed_columns = "-c 1 -S 2 -E 3 -g 4"
    log:
        f"{VAR_DIR}/01_Variants/Logs/{{sample}}.log",
    benchmark:
        f"{VAR_DIR}/01_Variants/Logs/{{sample}}.benchmark"
    threads: 8
    resources:
        mem_mb = 48000
    wrapper:
        f"{WRAPPER}/bio/vardict"
