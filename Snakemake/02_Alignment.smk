"""
PIPELINE: Alignment of DNA-seq data (no UMIs)
DESCRIPTION:
Read alignment with BWA-MEM2 (coordinate-sorted), duplicate marking with
Picard and base quality score recalibration with GATK. The final BAM is
indexed in 02_Alignment_QC.smk, shared with the UMI workflow.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Reads alignment. No -M: split alignments must stay supplementary (0x800),
# not secondary, so that Manta/Lumpy can use them as split-read evidence.
rule bwa_align:
    input:
        reads = [f"{QC_DIR}/{{sample}}_trimmed_R1.fastq.gz", f"{QC_DIR}/{{sample}}_trimmed_R2.fastq.gz"],
        idx = BWA_IDX
    output:
        temp(f"{ALIGN_DIR}/01_RawBAM/{{sample}}.bam")
    params:
        extra = r"-Y -R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
        sort = "samtools",
        sort_order = "coordinate"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_bwa.log"
    benchmark:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_bwa.benchmark"
    threads: 24
    resources:
        mem_mb = 64000
    wrapper:
        f"{WRAPPER}/bio/bwa-mem2/mem"

# Mark duplicates
rule mark_duplicates:
    input:
        bams = f"{ALIGN_DIR}/01_RawBAM/{{sample}}.bam"
    output:
        bam = temp(f"{ALIGN_DIR}/02_DedupBAM/{{sample}}_dedup.bam"),
        idx = temp(f"{ALIGN_DIR}/02_DedupBAM/{{sample}}_dedup.bam.bai"),
        metrics = f"{ALIGN_QC_DIR}/{{sample}}_dedup_stats.txt"
    params:
        extra = "--VALIDATION_STRINGENCY SILENT"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_markdup.log"
    resources:
        mem_mb = 36000
    wrapper:
        f"{WRAPPER}/bio/picard/markduplicates"

# Generate the recalibration table (restricted to the panel targets)
rule base_recalibrator:
    input:
        bam = f"{ALIGN_DIR}/02_DedupBAM/{{sample}}_dedup.bam",
        bai = f"{ALIGN_DIR}/02_DedupBAM/{{sample}}_dedup.bam.bai",
        ref = GENOME,
        dict = GENOME_DICT,
        fai = GENOME_FAI,
        known = [REFS['dbsnp'], REFS['mills'], REFS['high_conf']],
        bed = INTERVALS
    output:
        recal_table = f"{ALIGN_DIR}/{{sample}}_recal_data.tsv"
    params:
        extra = lambda wildcards, input: f"--intervals {input.bed}"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_bqsr.log"
    resources:
        mem_mb = 16000
    wrapper:
        f"{WRAPPER}/bio/gatk/baserecalibrator"

# Apply recalibration
rule apply_bqsr:
    input:
        bam = f"{ALIGN_DIR}/02_DedupBAM/{{sample}}_dedup.bam",
        bai = f"{ALIGN_DIR}/02_DedupBAM/{{sample}}_dedup.bam.bai",
        ref = GENOME,
        dict = GENOME_DICT,
        recal_table = f"{ALIGN_DIR}/{{sample}}_recal_data.tsv"
    output:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam"
    params:
        # The index is created by samtools as <bam>.bai in 02_Alignment_QC.smk
        extra = "--create-output-bam-index false"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_applybqsr.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/gatk/applybqsr"

# Inputs for the alignment MultiQC report (used by 02_Alignment_QC.smk)
ALIGN_MODE_QC = expand(f"{ALIGN_QC_DIR}/{{sample}}_dedup_stats.txt", sample=SAMPLES)
