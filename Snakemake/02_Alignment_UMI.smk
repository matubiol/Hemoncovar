"""
PIPELINE: Alignment and UMI-aware consensus calling for DNA-seq data
DESCRIPTION:
This pipeline takes paired-end FASTQ files, aligns them to a reference genome,
groups reads by UMI (unique molecular identifier), calls consensus reads, and
filters them for quality. It supports both duplex and simplex UMI strategies.
Steps without a Snakemake wrapper (fgbio FastqToBam, CorrectUmis, ZipperBams,
CallDuplexConsensusReads and the streamed re-alignments) use pinned conda
environments with the same fgbio/bwa-mem2/samtools versions as the wrappers.
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
UMI_STRATEGY = config['umi']['strat']  # "adjacency" (simplex), "duplex", or "none"
# Optional: path to a file with one expected UMI sequence per line, one entry per UMI
UMI_WHITELIST = config['umi']['list'] if config['umi']['list'] else None
UMI_DIR = f"{ALIGN_DIR}"
RG = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Normalizes each UMI segment to a fixed length (if necessary)
rule normalize_umi_length:
    input:
        r1 = f"{QC_DIR}/{{sample}}_trimmed_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}_trimmed_R2.fastq.gz"
    output:
        r1 = temp(f"{QC_DIR}/{{sample}}_umi_norm_R1.fastq.gz"),
        r2 = temp(f"{QC_DIR}/{{sample}}_umi_norm_R2.fastq.gz")
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_normalize_umi.log"
    threads: 2
    resources:
        mem_mb = 2000
    conda:
        ENVS['python']
    script:
        PRE['normalize_umi']

# Convert to unmapped BAM, extracting the UMI from the read name into a proper SAM tag (RX).
if UMI_WHITELIST:
    # Variant used when the UMI set is known and limited.
    rule fastq_to_ubam:
        input:
            r1 = f"{QC_DIR}/{{sample}}_umi_norm_R1.fastq.gz",
            r2 = f"{QC_DIR}/{{sample}}_umi_norm_R2.fastq.gz",
            umi_list = UMI_WHITELIST
        output:
            bam = temp(f"{UMI_DIR}/00_UnmappedBAM/{{sample}}_unmapped.bam"),
            metrics = f"{ALIGN_QC_DIR}/{{sample}}_umi_correction_metrics.txt"
        log:
            f"{RES_DIR}/02_Alignment/Logs/{{sample}}_fastq_to_ubam.log"
        conda:
            ENVS['fgbio']
        threads: 4
        resources:
            mem_mb = 10000
        shell:
            """
            set -o pipefail
            fgbio -Xmx4g --compression 0 FastqToBam \
                --input {input.r1} {input.r2} \
                --extract-umis-from-read-names \
                --sample {wildcards.sample} \
                --library {wildcards.sample} \
                --output /dev/stdout \
                2> {log} \
            | fgbio -Xmx4g --compression 1 --async-io CorrectUmis \
                --input /dev/stdin \
                --output {output.bam} \
                --metrics {output.metrics} \
                --max-mismatches 1 \
                --min-distance 2 \
                --umi-files {input.umi_list} \
                --dont-store-original-umis \
                2>> {log}
            """
else:
    # Variant used when the UMI is (or is treated as) fully degenerate/random
    rule fastq_to_ubam:
        input:
            r1 = f"{QC_DIR}/{{sample}}_umi_norm_R1.fastq.gz",
            r2 = f"{QC_DIR}/{{sample}}_umi_norm_R2.fastq.gz"
        output:
            temp(f"{UMI_DIR}/00_UnmappedBAM/{{sample}}_unmapped.bam")
        log:
            f"{RES_DIR}/02_Alignment/Logs/{{sample}}_fastq_to_ubam.log"
        conda:
            ENVS['fgbio']
        threads: 4
        resources:
            mem_mb = 10000
        shell:
            """
            fgbio -Xmx8g FastqToBam \
                --input {input.r1} {input.r2} \
                --extract-umis-from-read-names \
                --sample {wildcards.sample} \
                --library {wildcards.sample} \
                --output {output} \
                2> {log}
            """

# Align the unmapped BAM, keeping the UMI tags (streamed: samtools fastq | bwa-mem2 | ZipperBams).
# No -M: split alignments must stay supplementary so that SV callers can use them.
# ZipperBams streams its input, so a small heap is enough and leaves memory to bwa-mem2.
rule align_ubam:
    input:
        ref = GENOME,
        idx = BWA_IDX,
        dict = GENOME_DICT,
        ubam = f"{UMI_DIR}/00_UnmappedBAM/{{sample}}_unmapped.bam"
    output:
        temp(f"{UMI_DIR}/01_MappedTagged/{{sample}}_mapped_tagged.bam")
    params:
        rg = RG
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_align.log"
    threads: 24
    conda:
        ENVS['umi_align']
    resources:
        mem_mb = 64000
    shell:
        """
        set -o pipefail
        samtools fastq {input.ubam} 2> {log} \
        | bwa-mem2 mem -t {threads} -p -K 150000000 -Y \
            -R '{params.rg}' \
            {input.ref} - 2>> {log} \
        | fgbio -Xmx4g ZipperBams \
            --input /dev/stdin \
            --unmapped {input.ubam} \
            --ref {input.ref} \
            --output {output} \
            2>> {log}
        """

# Sort by TEMPLATE-COORDINATE (required by GroupReadsByUmi)
rule sort_template_coordinate:
    input:
        f"{UMI_DIR}/01_MappedTagged/{{sample}}_mapped_tagged.bam"
    output:
        temp(f"{UMI_DIR}/02_TemplateSorted/{{sample}}_template_sorted.bam")
    params:
        extra = "--template-coordinate"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_template_sort.log"
    threads: 8
    resources:
        mem_mb = 16000
    wrapper:
        f"{WRAPPER}/bio/samtools/sort"

# Group reads by UMI + mapping position
rule group_reads_by_umi:
    input:
        f"{UMI_DIR}/02_TemplateSorted/{{sample}}_template_sorted.bam"
    output:
        bam = temp(f"{UMI_DIR}/03_Grouped/{{sample}}_grouped.bam"),
        hist = f"{ALIGN_QC_DIR}/{{sample}}_umi_family_sizes.txt"
    params:
        extra = f"--strategy {'Paired' if UMI_STRATEGY == 'duplex' else 'Adjacency'} --edits 1"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_group_umi.log"
    threads: 8
    resources:
        mem_mb = 20000
    wrapper:
        f"{WRAPPER}/bio/fgbio/groupreadsbyumi"

if UMI_STRATEGY == "duplex":
    # Call DUPLEX consensus reads (no wrapper available)
    rule call_duplex_consensus:
        input:
            f"{UMI_DIR}/03_Grouped/{{sample}}_grouped.bam"
        output:
            temp(f"{UMI_DIR}/04_Consensus/{{sample}}_consensus_unmapped.bam")
        params:
            min_reads = "1 1 0"
        log:
            f"{RES_DIR}/02_Alignment/Logs/{{sample}}_call_consensus.log"
        conda:
            ENVS['fgbio']
        threads: 8
        resources:
            mem_mb = 24000
        shell:
            """
            fgbio -Xmx20g CallDuplexConsensusReads \
                --input {input} \
                --output {output} \
                --min-reads {params.min_reads} \
                --error-rate-pre-umi 45 \
                --error-rate-post-umi 30 \
                --threads {threads} \
                2> {log}
            """
else:
    # Call SIMPLEX consensus reads
    rule call_simplex_consensus:
        input:
            f"{UMI_DIR}/03_Grouped/{{sample}}_grouped.bam"
        output:
            temp(f"{UMI_DIR}/04_Consensus/{{sample}}_consensus_unmapped.bam")
        params:
            extra = "--min-reads 1 --max-reads 100 --error-rate-pre-umi 45 --error-rate-post-umi 40"
        log:
            f"{RES_DIR}/02_Alignment/Logs/{{sample}}_call_consensus.log"
        threads: 8
        resources:
            mem_mb = 24000
        wrapper:
            f"{WRAPPER}/bio/fgbio/callmolecularconsensusreads"

# Re-align the consensus reads and reattach consensus tags (streamed, no wrapper)
rule realign_consensus:
    input:
        ref = GENOME,
        idx = BWA_IDX,
        dict = GENOME_DICT,
        consensus_ubam = f"{UMI_DIR}/04_Consensus/{{sample}}_consensus_unmapped.bam"
    output:
        temp(f"{UMI_DIR}/06_Zippered/{{sample}}_zippered.bam")
    params:
        rg = RG
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_realign_consensus.log"
    threads: 24
    conda:
        ENVS['umi_align']
    resources:
        mem_mb = 36000
    shell:
        """
        set -o pipefail
        samtools fastq {input.consensus_ubam} 2> {log} \
        | bwa-mem2 mem -t {threads} -p -K 150000000 -Y \
            -R '{params.rg}' \
            {input.ref} - 2>> {log} \
        | fgbio -Xmx4g ZipperBams \
            --input /dev/stdin \
            --unmapped {input.consensus_ubam} \
            --ref {input.ref} \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
            --output {output} \
            2>> {log}
        """

# Query-name sort before filtering (all reads of a template together)
rule sort_consensus_queryname:
    input:
        f"{UMI_DIR}/06_Zippered/{{sample}}_zippered.bam"
    output:
        temp(f"{UMI_DIR}/07_QuerySorted/{{sample}}_zippered_qsorted.bam")
    params:
        extra = "-n"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_qsort.log"
    threads: 8
    resources:
        mem_mb = 16000
    wrapper:
        f"{WRAPPER}/bio/samtools/sort"

# Filter consensus reads by quality; output is coordinate-sorted and indexed in 02_Alignment_QC.smk
rule filter_consensus_reads:
    input:
        f"{UMI_DIR}/07_QuerySorted/{{sample}}_zippered_qsorted.bam",
        ref = GENOME
    output:
        f"{ALIGN_DIR}/{{sample}}_final.bam"
    params:
        ref = str(GENOME),
        min_reads = [1, 1, 0] if UMI_STRATEGY == "duplex" else [1],
        min_base_quality = 15,
        extra = "--max-read-error-rate 0.05 --max-base-error-rate 0.1 --max-no-call-fraction 0.2 --sort-order Coordinate"
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}_filter_consensus.log"
    threads: 8
    resources:
        mem_mb = 16000
    wrapper:
        f"{WRAPPER}/bio/fgbio/filterconsensusreads"

# Inputs for the alignment MultiQC report (used by 02_Alignment_QC.smk)
ALIGN_MODE_QC = expand(f"{ALIGN_QC_DIR}/{{sample}}_umi_family_sizes.txt", sample=SAMPLES)
