"""
PIPELINE: Final BAM indexing and alignment QC (shared by UMI and non-UMI modes)
DESCRIPTION:
Indexes the final BAM (<sample>_final.bam.bai), collects alignment summary
metrics with Picard and aggregates the alignment QC with MultiQC. Mode-specific
QC files (duplication metrics or UMI family sizes) are provided by the
alignment module through ALIGN_MODE_QC.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
rule index_final_bam:
    input:
        f"{ALIGN_DIR}/{{sample}}_final.bam"
    output:
        f"{ALIGN_DIR}/{{sample}}_final.bam.bai"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_index.log"
    threads: 4
    wrapper:
        f"{WRAPPER}/bio/samtools/index"

# Alignment summary metrics
rule collect_metrics:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        ref = GENOME
    output:
        f"{ALIGN_QC_DIR}/{{sample}}_alignment_metrics.txt"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_alignment_metrics.log"
    resources:
        mem_mb = 4000
    wrapper:
        f"{WRAPPER}/bio/picard/collectalignmentsummarymetrics"

# Aggregate QC results with MultiQC
rule multiqc_align:
    input:
        expand(f"{ALIGN_QC_DIR}/{{sample}}_alignment_metrics.txt", sample=SAMPLES),
        ALIGN_MODE_QC
    output:
        f"{ALIGN_QC_DIR}/multiqc_report.html"
    params:
        extra = "--force",
        use_input_files_only = True
    log:
        f"{ALIGN_QC_DIR}/multiqc.log"
    wrapper:
        f"{WRAPPER}/bio/multiqc"
