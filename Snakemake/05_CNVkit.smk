"""
PIPELINE: CNVkit Copy Number Variation Analysis
DESCRIPTION:
Step-wise CNVkit workflow (coverage -> reference -> fix -> segment -> call).
The 'batch' wrapper is not used: it recomputes the coverage of every normal
each time a reference is built, which makes the leave-one-out PoN QC
quadratic. Coverage is computed once per sample and reused.

  Panel of normals (build mode, written to PON_DIR/cnvkit):
    access.bed, targets.bed, antitargets.bed, reference.cnn
  Per tumour sample:
    coverage -> fix (vs reference) -> segment -> call (per-sample purity/sex,
    B-allele frequencies from heterozygous SNPs genotyped with bcftools)
    -> scatter / diagram -> gene/cytoband annotation
  PoN QC (build mode, pon_qc): every normal is fixed, segmented and called
    against a reference built without it; sex is inferred with 'cnvkit.py sex'.

Wrappers: cnvkit target/antitarget/call/diagram, bcftools mpileup/call.
Pinned conda environment (cnvkit 0.9.14, as the wrappers): access, coverage,
reference, fix, segment, scatter, sex.
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
CNVKIT_PON = f"{PON_DIR}/cnvkit"

def sample_sex_opt(wildcards):
    """'--sample-sex <sex>' for tumours or normals, when the sex is known."""
    sample = getattr(wildcards, "sample", None) or wildcards.normal
    sex = SAMPLE_INFO[sample]["sex"]
    return f"--sample-sex {sex}" if sex else ""

def other_normals(wildcards, kind):
    """Coverage files of all normals except the one being evaluated (leave-one-out)."""
    return [f"{CNVkit_DIR}/01_Coverage/{n}.{kind}coverage.cnn" for n in NORMALS if n != wildcards.normal]

# =============================================================================
# PANEL OF NORMALS (build mode)
# =============================================================================
if PON_MODE == "build":
    # Accessible regions of the reference (excludes N stretches) - no wrapper
    rule cnvkit_access:
        input:
            fasta = GENOME
        output:
            f"{CNVKIT_PON}/access.bed"
        log:
            f"{CNVKIT_PON}/Logs/access.log"
        conda:
            ENVS['cnvkit']
        shell:
            "cnvkit.py access {input.fasta} -o {output} 2> {log}"

    # Target bins (split into ~uniform bins, annotated with gene names)
    rule cnvkit_target:
        input:
            bed = INTERVALS,
            annotate = REFS['gtf']
        output:
            bed = f"{CNVKIT_PON}/targets.bed"
        params:
            extra = "--split"
        log:
            f"{CNVKIT_PON}/Logs/target.log"
        wrapper:
            f"{WRAPPER}/bio/cnvkit/target"

    # Off-target (antitarget) bins: genome-wide backbone for CNV
    rule cnvkit_antitarget:
        input:
            bed = f"{CNVKIT_PON}/targets.bed",
            accessible = f"{CNVKIT_PON}/access.bed"
        output:
            bed = f"{CNVKIT_PON}/antitargets.bed"
        log:
            f"{CNVKIT_PON}/Logs/antitarget.log"
        wrapper:
            f"{WRAPPER}/bio/cnvkit/antitarget"

    # Reference from the normal samples (sex of each normal inferred by CNVkit) - no wrapper
    rule cnvkit_reference:
        input:
            cnn = expand(f"{CNVkit_DIR}/01_Coverage/{{normal}}.{{kind}}coverage.cnn",
                         normal=NORMALS, kind=["target", "antitarget"]),
            fasta = GENOME,
            fai = GENOME_FAI
        output:
            f"{CNVKIT_PON}/reference.cnn"
        log:
            f"{CNVKIT_PON}/Logs/reference.log"
        conda:
            ENVS['cnvkit']
        resources:
            mem_mb = 8000
        shell:
            "cnvkit.py reference {input.cnn} -f {input.fasta} -o {output} 2> {log}"

# =============================================================================
# PER-SAMPLE RULES (normals and tumours)
# =============================================================================
# Target and antitarget coverage - no wrapper
rule cnvkit_coverage:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_final.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        targets = f"{CNVKIT_PON}/targets.bed",
        antitargets = f"{CNVKIT_PON}/antitargets.bed",
        fasta = GENOME
    output:
        target = f"{CNVkit_DIR}/01_Coverage/{{sample}}.targetcoverage.cnn",
        antitarget = f"{CNVkit_DIR}/01_Coverage/{{sample}}.antitargetcoverage.cnn"
    log:
        f"{CNVkit_DIR}/01_Coverage/Logs/{{sample}}.log"
    conda:
        ENVS['cnvkit']
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        """
        cnvkit.py coverage {input.bam} {input.targets} -f {input.fasta} -p {threads} -o {output.target} 2> {log}
        cnvkit.py coverage {input.bam} {input.antitargets} -f {input.fasta} -p {threads} -o {output.antitarget} 2>> {log}
        """

# Heterozygous SNP genotyping at common SNPs within the targets, for B-allele frequencies
rule het_snps_mpileup:
    input:
        alignments = f"{ALIGN_DIR}/{{sample}}_final.bam",
        index = f"{ALIGN_DIR}/{{sample}}_final.bam.bai",
        ref = GENOME,
        regions = INTERVALS,
        targets = REFS['common_snps']
    output:
        temp(f"{CNVkit_DIR}/00_HetSNPs/{{sample}}.pileup.bcf")
    params:
        extra = "--annotate FORMAT/AD,FORMAT/DP --max-depth 10000 --min-MQ 20 --min-BQ 20"
    log:
        f"{CNVkit_DIR}/00_HetSNPs/Logs/{{sample}}_mpileup.log"
    wrapper:
        f"{WRAPPER}/bio/bcftools/mpileup"

rule het_snps_call:
    input:
        f"{CNVkit_DIR}/00_HetSNPs/{{sample}}.pileup.bcf"
    output:
        f"{CNVkit_DIR}/00_HetSNPs/{{sample}}.vcf.gz"
    params:
        caller = "-m",
        extra = "--variants-only --skip-variants indels"
    log:
        f"{CNVkit_DIR}/00_HetSNPs/Logs/{{sample}}_call.log"
    wrapper:
        f"{WRAPPER}/bio/bcftools/call"

# =============================================================================
# TUMOUR SAMPLES
# =============================================================================
# Normalize against the reference - no wrapper
rule cnvkit_fix:
    input:
        target = f"{CNVkit_DIR}/01_Coverage/{{sample}}.targetcoverage.cnn",
        antitarget = f"{CNVkit_DIR}/01_Coverage/{{sample}}.antitargetcoverage.cnn",
        reference = f"{CNVKIT_PON}/reference.cnn",
        pon_check = PON_CHECK
    output:
        f"{CNVkit_DIR}/02_Fix/{{sample}}.cnr"
    log:
        f"{CNVkit_DIR}/02_Fix/Logs/{{sample}}.log"
    conda:
        ENVS['cnvkit']
    resources:
        mem_mb = 8000
    shell:
        """
        cnvkit.py fix {input.target} {input.antitarget} -r {input.reference} \
            -i {wildcards.sample} -o {output} 2> {log}
        """

# Segmentation (CBS) - no wrapper
rule cnvkit_segment:
    input:
        f"{CNVkit_DIR}/02_Fix/{{sample}}.cnr"
    output:
        f"{CNVkit_DIR}/03_Segment/{{sample}}.cns"
    log:
        f"{CNVkit_DIR}/03_Segment/Logs/{{sample}}.log"
    conda:
        ENVS['cnvkit']
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        """
        cnvkit.py segment {input} --drop-low-coverage -p {threads} -o {output} 2> {log}
        """

# Absolute copy number with per-sample purity, sex and B-allele frequencies.
# '--filter' goes through 'extra': the 'filter' param of the wrapper is not
# formatted correctly (it would pass the literal string '{filter}').
rule cnvkit_call:
    input:
        segment = f"{CNVkit_DIR}/03_Segment/{{sample}}.cns",
        vcf = f"{CNVkit_DIR}/00_HetSNPs/{{sample}}.vcf.gz"
    output:
        segment = f"{CNVkit_DIR}/04_Call/{{sample}}.call.cns"
    params:
        # Without a purity value no purity rescaling is applied
        purity = lambda wildcards: SAMPLE_INFO[wildcards.sample]["purity"],
        ploidy = 2,
        extra = lambda wildcards: (f"--filter ampdel --drop-low-coverage --sample-id {wildcards.sample} "
                                   f"{sample_sex_opt(wildcards)}")
    log:
        f"{CNVkit_DIR}/04_Call/Logs/{{sample}}_call.log"
    wrapper:
        f"{WRAPPER}/bio/cnvkit/call"

# Scatter plot with B-allele frequencies - no wrapper
rule cnvkit_scatter:
    input:
        cnr = f"{CNVkit_DIR}/02_Fix/{{sample}}.cnr",
        cns = f"{CNVkit_DIR}/04_Call/{{sample}}.call.cns",
        vcf = f"{CNVkit_DIR}/00_HetSNPs/{{sample}}.vcf.gz"
    output:
        f"{CNVkit_DIR}/04_Call/{{sample}}_scatter.pdf"
    log:
        f"{CNVkit_DIR}/04_Call/Logs/{{sample}}_scatter.log"
    conda:
        ENVS['cnvkit']
    shell:
        """
        cnvkit.py scatter {input.cnr} -s {input.cns} -v {input.vcf} -i {wildcards.sample} \
            -o {output} 2> {log}
        """

# Chromosome diagram
rule cnvkit_diagram:
    input:
        f"{CNVkit_DIR}/02_Fix/{{sample}}.cnr",
        f"{CNVkit_DIR}/04_Call/{{sample}}.call.cns"
    output:
        f"{CNVkit_DIR}/04_Call/{{sample}}_diagram.pdf"
    params:
        extra = sample_sex_opt
    log:
        f"{CNVkit_DIR}/04_Call/Logs/{{sample}}_diagram.log"
    wrapper:
        f"{WRAPPER}/bio/cnvkit/diagram"

# Gene and cytoband annotation of all segments
rule annotate_coord_cnvkit:
    input:
        call = f"{CNVkit_DIR}/04_Call/{{sample}}.call.cns",
        cyto = REFS['cyto']
    params:
        format = "cnvkit",
        min_probes = 0,
        sample = "{sample}",
        sex = lambda wildcards: SAMPLE_INFO[wildcards.sample]["sex"]
    output:
        tsv = f"{CNVkit_DIR}/05_Annotated/{{sample}}.ann.tsv"
    log:
        f"{CNVkit_DIR}/05_Annotated/Logs/{{sample}}.log"
    conda:
        ENVS['python']
    script:
        SRC['CNV_gene_ann']

# =============================================================================
# PoN QC: leave-one-out evaluation of the normals (build mode)
# =============================================================================
if PON_QC:
    rule cnvkit_loo_reference:
        input:
            target = lambda wildcards: other_normals(wildcards, "target"),
            antitarget = lambda wildcards: other_normals(wildcards, "antitarget"),
            fasta = GENOME
        output:
            temp(f"{CNVkit_DIR}/PoN_QC/{{normal}}.loo_reference.cnn")
        log:
            f"{CNVkit_DIR}/PoN_QC/Logs/{{normal}}_reference.log"
        conda:
            ENVS['cnvkit']
        resources:
            mem_mb = 8000
        shell:
            "cnvkit.py reference {input.target} {input.antitarget} -f {input.fasta} -o {output} 2> {log}"

    rule cnvkit_loo_call:
        input:
            target = f"{CNVkit_DIR}/01_Coverage/{{normal}}.targetcoverage.cnn",
            antitarget = f"{CNVkit_DIR}/01_Coverage/{{normal}}.antitargetcoverage.cnn",
            reference = f"{CNVkit_DIR}/PoN_QC/{{normal}}.loo_reference.cnn"
        output:
            cnr = f"{CNVkit_DIR}/PoN_QC/{{normal}}.cnr",
            cns = temp(f"{CNVkit_DIR}/PoN_QC/{{normal}}.cns"),
            call = f"{CNVkit_DIR}/PoN_QC/{{normal}}.call.cns"
        params:
            sex = sample_sex_opt
        log:
            f"{CNVkit_DIR}/PoN_QC/Logs/{{normal}}_call.log"
        conda:
            ENVS['cnvkit']
        threads: 4
        resources:
            mem_mb = 8000
        shell:
            """
            cnvkit.py fix {input.target} {input.antitarget} -r {input.reference} \
                -i {wildcards.normal} -o {output.cnr} 2> {log}
            cnvkit.py segment {output.cnr} --drop-low-coverage -p {threads} -o {output.cns} 2>> {log}
            cnvkit.py call {output.cns} --filter ampdel --drop-low-coverage {params.sex} \
                -o {output.call} 2>> {log}
            """

    # Inferred sex and X/Y log2 ratios of every normal
    rule cnvkit_loo_sex:
        input:
            expand(f"{CNVkit_DIR}/PoN_QC/{{normal}}.cnr", normal=NORMALS)
        output:
            f"{CNVkit_DIR}/PoN_QC/cnvkit_sex.tsv"
        log:
            f"{CNVkit_DIR}/PoN_QC/Logs/sex.log"
        conda:
            ENVS['cnvkit']
        shell:
            "cnvkit.py sex {input} -o {output} 2> {log}"
