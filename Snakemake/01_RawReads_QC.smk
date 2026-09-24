"""
PIPELINE: Raw Reads Quality Control
DESCRIPTION:
This Snakemake workflow performs quality control on raw DNA-seq reads,
including FastQC for quality assessment, Fastp for trimming and filtering,
and MultiQC for aggregating QC reports. Samples sequenced across several
lanes are concatenated first, so no lane is silently dropped.
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
MERGED_DIR = f"{RES_DIR}/01_RawReads_QC/00_MergedLanes"

def get_raw(wildcards, read=None):
    """Return the FASTQ to use for a sample/read: the file itself if the sample
    has a single lane, otherwise the lane-merged FASTQ (built by merge_lanes)."""
    read = read or wildcards.read
    files = SAMPLE_INFO[wildcards.sample][read]
    if len(files) == 1:
        return files[0]
    return f"{MERGED_DIR}/{wildcards.sample}_{read}.fastq.gz"

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Concatenate lanes (gzip streams can be concatenated directly)
rule merge_lanes:
    input:
        lambda wildcards: SAMPLE_INFO[wildcards.sample][wildcards.read]
    output:
        temp(f"{MERGED_DIR}/{{sample}}_{{read}}.fastq.gz")
    shell:
        "cat {input} > {output}"

# FastQC for raw reads
rule fastqc:
    input:
        get_raw
    output:
        html = f"{RES_DIR}/01_RawReads_QC/01_FastQC/{{sample}}_{{read}}_fastqc.html",
        zip  = f"{RES_DIR}/01_RawReads_QC/01_FastQC/{{sample}}_{{read}}_fastqc.zip"
    log:
        f"{RES_DIR}/01_RawReads_QC/01_FastQC/Logs/{{sample}}_{{read}}.log"
    threads: 4
    resources:
        mem_mb = 4096
    wrapper:
        f"{WRAPPER}/bio/fastqc"

# Fastp for adapter trimming and quality filtering
rule fastp:
    input:
        sample = lambda wildcards: [get_raw(wildcards, "R1"), get_raw(wildcards, "R2")]
    output:
        trimmed = [f"{QC_DIR}/{{sample}}_trimmed_R1.fastq.gz", f"{QC_DIR}/{{sample}}_trimmed_R2.fastq.gz"],
        html = f"{QC_DIR}/{{sample}}_fastp_report.html",
        json = f"{QC_DIR}/{{sample}}_fastp_report.json"
    params:
        extra = "--detect_adapter_for_pe --length_required 35 --cut_right"
    log:
        f"{QC_DIR}/Logs/{{sample}}.log"
    threads: 8
    wrapper:
        f"{WRAPPER}/bio/fastp"

# MultiQC to aggregate FastQC and Fastp reports
rule raw_reads_multiqc:
    input:
        expand(f"{RES_DIR}/01_RawReads_QC/01_FastQC/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        expand(f"{QC_DIR}/{{sample}}_fastp_report.json", sample=SAMPLES)
    output:
        f"{RES_DIR}/01_RawReads_QC/03_MultiQC/MultiQC_report.html"
    params:
        extra = "--force",
        use_input_files_only = True
    log:
        f"{RES_DIR}/01_RawReads_QC/03_MultiQC/multiqc.log"
    wrapper:
        f"{WRAPPER}/bio/multiqc"
