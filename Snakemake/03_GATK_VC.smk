"""
PIPELINE: GATK Small Variation Analysis
DESCRIPTION:
Somatic variant calling with GATK Mutect2 (tumour-only) followed by the
GATK filtering chain: read-orientation model, contamination estimate and
FilterMutectCalls. All steps use Snakemake wrappers (GATK 4.6.2.0). The filtered VCF (all records, FILTER column filled) is
then normalized, annotated and reported by 03_SmallVariants_Common.smk.
Mutect2 sensitivity is configurable through 'mutect_mode'.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Dynamically resolve Mutect2 execution parameters based on configuration modes.
def get_mutect_params(wildcards):
    mode = config.get("mutect_mode", "moderate")
    common = "--pair-hmm-implementation FASTEST_AVAILABLE --genotype-pon-sites true"

    # Define extra flags as a clean continuous string (no commas at line ends)
    extra = (
        "--use-pdhmm true --use-pdhmm-overlap-optimization true "
        "--disable-read-filter SoftClippedReadFilter "
        "--disable-adaptive-pruning true --max-mnp-distance 10"
    )

    if mode == "default":
        return common

    elif mode == "moderate":
        moderated_flags = [
            common,
            extra,
            "--initial-tumor-lod 0.5",
            "--tumor-lod-to-emit 0.5",
            "--kmer-size 10",
            "--kmer-size 25"
        ]
        # Join all elements into a single clean space-separated string
        return " ".join(moderated_flags)

    elif mode == "relaxed":
        relaxed_flags = [
            common,
            extra,
            "--initial-tumor-lod 0.05",
            "--tumor-lod-to-emit 0.05",
            "--kmer-size 10",
            "--kmer-size 15",
            "--kmer-size 25",
            "--active-probability-threshold 0.002"
        ]
        return " ".join(relaxed_flags)

    return common

# Run Mutect2 for somatic variant calling
rule mutect2:
    input:
        map = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        fasta = GENOME,
        dict = GENOME_DICT,
        germline = REFS['gnomad'],
        pon = REFS['pon'],
        intervals = INTERVALS
    output:
        vcf = temp(f"{GATK_DIR}/01_Variants/{{sample}}.unfiltered.vcf.gz"),
        tbi = temp(f"{GATK_DIR}/01_Variants/{{sample}}.unfiltered.vcf.gz.tbi"),
        stats = temp(f"{GATK_DIR}/01_Variants/{{sample}}.unfiltered.vcf.gz.stats"),
        f1r2 = temp(f"{GATK_DIR}/01_Variants/{{sample}}_f1r2.tar.gz")
    params:
        extra = get_mutect_params
    log:
        f"{GATK_DIR}/01_Variants/Logs/{{sample}}.log",
    benchmark:
        f"{GATK_DIR}/01_Variants/Logs/{{sample}}.benchmark"
    threads: 8
    resources:
        mem_mb = 48000
    wrapper:
        f"{WRAPPER}/bio/gatk/mutect"

# Learn the read-orientation (FFPE / oxidation) artefact model
rule learn_read_orientation:
    input:
        f1r2 = f"{GATK_DIR}/01_Variants/{{sample}}_f1r2.tar.gz"
    output:
        temp(f"{GATK_DIR}/01_Variants/{{sample}}_read-orientation-model.tar.gz")
    log:
        f"{GATK_DIR}/01_Variants/Logs/{{sample}}_orientation.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/learnreadorientationmodel"

# Pileups at common SNPs within the panel, for the contamination estimate
rule get_pileup_summaries:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        intervals = REFS['common_snps'],
        variants = REFS['common_snps'],
        bed = INTERVALS
    output:
        temp(f"{GATK_DIR}/01_Variants/{{sample}}_pileups.table")
    params:
        # Restrict the common SNPs to the panel targets
        extra = lambda wildcards, input: f"--intervals {input.bed} --interval-set-rule INTERSECTION"
    log:
        f"{GATK_DIR}/01_Variants/Logs/{{sample}}_pileups.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/getpileupsummaries"

# Estimate cross-sample contamination.
# Note: the wrapper does not expose --tumor-segmentation; FilterMutectCalls then
# works without the tumour segmentation (minor impact on small panels).
rule calculate_contamination:
    input:
        tumor = f"{GATK_DIR}/01_Variants/{{sample}}_pileups.table"
    output:
        temp(f"{GATK_DIR}/01_Variants/{{sample}}_contamination.table")
    log:
        f"{GATK_DIR}/01_Variants/Logs/{{sample}}_contamination.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/calculatecontamination"

# Apply the Mutect2 filters (germline, PoN, strand/orientation bias, contamination...)
rule filter_mutect_calls:
    input:
        vcf = f"{GATK_DIR}/01_Variants/{{sample}}.unfiltered.vcf.gz",
        tbi = f"{GATK_DIR}/01_Variants/{{sample}}.unfiltered.vcf.gz.tbi",
        stats = f"{GATK_DIR}/01_Variants/{{sample}}.unfiltered.vcf.gz.stats",
        ref = GENOME,
        dict = GENOME_DICT,
        f1r2 = f"{GATK_DIR}/01_Variants/{{sample}}_read-orientation-model.tar.gz",
        contamination_table = f"{GATK_DIR}/01_Variants/{{sample}}_contamination.table"
    output:
        vcf = temp(f"{GATK_DIR}/01_Variants/{{sample}}.filtered.vcf.gz"),
        tbi = temp(f"{GATK_DIR}/01_Variants/{{sample}}.filtered.vcf.gz.tbi"),
        stats = f"{GATK_DIR}/01_Variants/{{sample}}.filtering_stats.tsv"
    params:
        extra = lambda wildcards, input, output: f"--stats {input.stats} --filtering-stats {output.stats}"
    log:
        f"{GATK_DIR}/01_Variants/Logs/{{sample}}_filter.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/filtermutectcalls"
