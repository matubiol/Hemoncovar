"""
PIPELINE: Small variant normalization, annotation and report (shared by all callers)
DESCRIPTION:
For every small variant caller ({caller} = GATK | VarDict):
  1. Normalizes the caller VCF with bcftools wrappers (contig header from the
     .fai, left-alignment, multiallelic split, sorting, bgzip + tabix).
  2. Annotates it with SnpEff/SnpSift wrappers (functional effect, variant
     type, ClinVar) and extracts the report fields (SnpSift extractFields).
  3. Builds the per-sample TSV report (03_SnpEff_report.py).
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
# Caller-specific VCF that enters the normalization step
CALLER_VCF = {
    "GATK": f"{GATK_DIR}/01_Variants/{{sample}}.filtered.vcf.gz",
    "VarDict": f"{VAR_DIR}/01_Variants/{{sample}}.raw.vcf",
}

# Fields extracted from the annotated VCF for the report
REPORT_FIELDS = [
    "CHROM", "POS", "REF", "ALT", "FILTER",
    "ANN[0].GENE", "ANN[0].FEATUREID", "ANN[0].RANK", "ANN[0].HGVS_C", "ANN[0].HGVS_P",
    "ANN[0].EFFECT", "ANN[0].IMPACT",
    "GEN[0].DP", "GEN[0].AD[1]", "GEN[0].AF",
    "CLNSIG", "VARTYPE",
]

# =============================================================================
# PIPELINE RULES
# =============================================================================
# 1a. Add/refresh contig lines from the .fai (VarDict VCFs have none)
rule reheader_small_vcf:
    input:
        lambda wildcards: CALLER_VCF[wildcards.caller].format(sample=wildcards.sample),
        fai = GENOME_FAI
    output:
        temp(f"{SMALL_DIR}/{{caller}}/01_Variants/{{sample}}.reheader.bcf")
    params:
        extra = lambda wildcards, input: f"--fai {input.fai}"
    log:
        f"{SMALL_DIR}/{{caller}}/01_Variants/Logs/{{sample}}_reheader.log"
    wrapper:
        f"{WRAPPER}/bio/bcftools/reheader"

# 1b. Left-align and split multiallelic records
rule norm_small_vcf:
    input:
        f"{SMALL_DIR}/{{caller}}/01_Variants/{{sample}}.reheader.bcf",
        ref = GENOME,
        fai = GENOME_FAI
    output:
        temp(f"{SMALL_DIR}/{{caller}}/01_Variants/{{sample}}.norm.bcf")
    params:
        extra = "--multiallelics -any --check-ref w"
    log:
        f"{SMALL_DIR}/{{caller}}/01_Variants/Logs/{{sample}}_norm.log"
    wrapper:
        f"{WRAPPER}/bio/bcftools/norm"

# 1c. Sort, compress and index
rule sort_small_vcf:
    input:
        f"{SMALL_DIR}/{{caller}}/01_Variants/{{sample}}.norm.bcf"
    output:
        f"{SMALL_DIR}/{{caller}}/01_Variants/{{sample}}.vcf.gz",
        f"{SMALL_DIR}/{{caller}}/01_Variants/{{sample}}.vcf.gz.tbi"
    params:
        extra = "--write-index=tbi"
    log:
        f"{SMALL_DIR}/{{caller}}/01_Variants/Logs/{{sample}}_sort.log"
    resources:
        mem_mb = 4000
    wrapper:
        f"{WRAPPER}/bio/bcftools/sort"

# 2a. SnpEff database (downloaded only if it is not already present in refs.snpeff_data)
rule snpeff_download:
    output:
        directory(SNPEFF_DB)
    params:
        reference = config['snpeff_db']
    log:
        f"{RES_DIR}/00_Reference/Logs/snpeff_download.log"
    resources:
        mem_mb = 4000
    wrapper:
        f"{WRAPPER}/bio/snpeff/download"

# 2b. Functional annotation (restricted to the panel targets)
rule snpeff_annotate:
    input:
        calls = f"{SMALL_DIR}/{{caller}}/01_Variants/{{sample}}.vcf.gz",
        db = SNPEFF_DB,
        bed = INTERVALS
    output:
        calls = temp(f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_ann.vcf")
    params:
        extra = lambda wildcards, input: f"-filterInterval {input.bed}"
    log:
        f"{SMALL_DIR}/{{caller}}/02_Annotated/Logs/{{sample}}_snpeff.log"
    resources:
        mem_mb = 18000
    wrapper:
        f"{WRAPPER}/bio/snpeff/annotate"

# 2c. Variant type (SNP, MNP, INS, DEL, MIXED)
rule snpsift_vartype:
    input:
        vcf = f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_ann.vcf"
    output:
        vcf = temp(f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_ann_vartype.vcf")
    log:
        f"{SMALL_DIR}/{{caller}}/02_Annotated/Logs/{{sample}}_vartype.log"
    resources:
        mem_mb = 4000
    wrapper:
        f"{WRAPPER}/bio/snpsift/varType"

# 2d. ClinVar significance
rule snpsift_clinvar:
    input:
        call = f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_ann_vartype.vcf",
        database = REFS['clinvar']
    output:
        call = temp(f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_ann_clinvar.vcf")
    params:
        extra = "-info CLNSIG"
    log:
        f"{SMALL_DIR}/{{caller}}/02_Annotated/Logs/{{sample}}_clinvar.log"
    resources:
        mem_mb = 8000
    wrapper:
        f"{WRAPPER}/bio/snpsift/annotate"

# 2e. Table for the report (SnpSift extractFields has no wrapper)
rule snpsift_extract_fields:
    input:
        f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_ann_clinvar.vcf"
    output:
        temp(f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_ann_clinvar.tsv")
    params:
        fields = " ".join(f'"{f}"' for f in REPORT_FIELDS)
    log:
        f"{SMALL_DIR}/{{caller}}/02_Annotated/Logs/{{sample}}_extract.log"
    conda:
        ENVS['snpsift']
    resources:
        mem_mb = 4000
    shell:
        """
        SnpSift extractFields -Xmx3g {input} {params.fields} > {output} 2> {log}
        """

# 3. Per-sample report
rule report_small:
    input:
        vars = f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_ann_clinvar.tsv"
    params:
        min_depth = config.get("vc_min_depth", 8),
        min_vaf = config.get("vc_min_vaf", 0.0),
        keep_filtered = config.get("report_keep_filtered", False)
    output:
        report = f"{SMALL_DIR}/{{caller}}/02_Annotated/{{sample}}_report.tsv"
    log:
        f"{SMALL_DIR}/{{caller}}/02_Annotated/Logs/{{sample}}_report.log"
    conda:
        ENVS['python']
    script:
        SRC['snpeff_report']
