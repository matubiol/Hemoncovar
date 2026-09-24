"""
PIPELINE: Hemoncovar - Haemato-oncology somatic variant, SV/fusion & CNV calling
DESCRIPTION:
    High-throughput pipeline for haemato-oncology NGS capture panels
    (designed around the HealthInCode Haematology OncoKit HS2: 89 genes,
    26 fusion genes with intronic baits, genome-wide CNV backbone, dual UMIs).
    It integrates raw data QC, alignment (optionally UMI-aware), small
    variant calling (GATK Mutect2, VarDict), structural variant / fusion
    detection (Manta, Lumpy) and copy number profiling (CNVkit, GATK CNV).
USAGE:
    snakemake -s Snakefile.smk --configfile config.yml --profile slurm \
        --software-deployment-method conda
VERSION:
    1.0.0
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
from pathlib import Path
import re
import pandas as pd

# Default config file (can be overridden with --configfile)
configfile: "config.yml"

# Set up paths based on config
PROJ_DIR = config['proj']['proj_dir']
RAW_DIR = config['proj']["raw_dir"]
RES_DIR = f"{PROJ_DIR}/results"
REF_DIR = Path(config['proj']['ref_dir'])
SRC_DIR = Path(config['proj']['src_dir'])

PIPELINE_VERSION = "1.0.0"

# Resolve absolute paths for references, modules and pre/post-processing scripts
REFS = {key: REF_DIR / path for key, path in config["refs"].items()}
SMK = {key: SRC_DIR / path for key, path in config["smk"].items()}
PRE = {key: SRC_DIR / path for key, path in config["pre"].items()}
SRC = {key: SRC_DIR / path for key, path in config["post"].items()}

# =============================================================================
# SOFTWARE
# =============================================================================
# Tools with a Snakemake wrapper use 'wrapper: f"{WRAPPER}/bio/..."' (pinned release).
# The remaining rules use the conda environments in envs/, pinned to the same
# tool versions as that wrapper release.
WRAPPER = config.get("wrapper_version", "v9.18.0")
ENVS = {p.stem: str(p) for p in (Path(workflow.basedir) / "envs").glob("*.yaml")}

# Reference companion files, declared as inputs so that missing indices fail early
GENOME = REFS['genome']
GENOME_FAI = f"{GENOME}.fai"
GENOME_DICT = str(GENOME.with_suffix(".dict"))
BWA_IDX = [f"{GENOME}{ext}" for ext in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")]
SNPEFF_DB = str(REFS['snpeff_data'] / config['snpeff_db'])

# =============================================================================
# INTERVALS CHROMOSOME NAMING CHECK
# =============================================================================
# Panel bed files and reference fastas often disagree on chromosome naming
CANONICAL_CHROM_RE = re.compile(r'^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$', re.IGNORECASE)

def fasta_has_chr_prefix(ref_fasta):
    """Scan fasta headers for the first canonical chromosome to check naming,
    ignoring scaffolds/contigs that don't match chr1-22/X/Y/M."""
    with open(ref_fasta) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if CANONICAL_CHROM_RE.match(name):
                    return name.lower().startswith("chr")
    raise ValueError(
        f"Could not find a canonical chromosome header (chr1-22/X/Y/M or 1-22/X/Y/M) "
        f"in reference fasta: {ref_fasta}"
    )

def bed_has_chr_prefix(bed_file):
    """Scan bed intervals for the first canonical chromosome to check naming,
    ignoring scaffolds/contigs that don't match chr1-22/X/Y/M."""
    with open(bed_file) as f:
        for line in f:
            if line.startswith(("track", "browser", "#")):
                continue
            chrom = line.split("\t")[0]
            if CANONICAL_CHROM_RE.match(chrom):
                return chrom.lower().startswith("chr")
    raise ValueError(
        f"Could not find a canonical chromosome interval (chr1-22/X/Y/M or 1-22/X/Y/M) "
        f"in bed file: {bed_file}"
    )

RAW_BED = config['proj']['intervals']
NEEDS_BED_HARMONIZATION = fasta_has_chr_prefix(REFS['genome']) != bed_has_chr_prefix(RAW_BED)

if NEEDS_BED_HARMONIZATION:
    print(f"[WARNING] Chromosome naming mismatch detected between '{RAW_BED}' and '{REFS['genome']}'. "
          f"The intervals bed will be harmonized in 00_RefPrep.smk.")

# =============================================================================
# SAMPLES
# =============================================================================
# Two ways of defining samples:
#   1. A sample sheet (config 'samplesheet'), recommended. Tab-separated with
#      columns: sample, fastq_1, fastq_2 [, type, sex, purity]
#      - fastq_1/fastq_2: one or several comma-separated FASTQs (one per lane),
#        absolute or relative to raw_dir. Lanes are concatenated in the order given.
#      - type:   'tumor' (default) or 'normal' (used for the CNV panels of normals)
#      - sex:    'male' / 'female' (optional, used for CNV calling)
#      - purity: tumour fraction in (0, 1] (optional, e.g. blast % / 100)
#   2. Auto-discovery in raw_dir (no sample sheet). All lanes of a sample are
#      kept. Sample names may contain underscores. Samples listed in
#      config 'pon_samples' are treated as normals.
ILLUMINA_FASTQ_RE = re.compile(r"^(?P<sample>.+?)_S\d+_L\d{3}_(?P<read>R[12])_001\.fastq\.gz$")
SIMPLE_FASTQ_RE = re.compile(r"^(?P<sample>.+?)_(?P<read>R[12])(?:_001)?\.fastq\.gz$")
VALID_SEX = {"male": "male", "m": "male", "female": "female", "f": "female", "": ""}

def _resolve_fastq(path):
    p = Path(path.strip())
    return str(p if p.is_absolute() else Path(RAW_DIR) / p)

def discover_samples(raw_dir):
    """Group FASTQs by sample and read, keeping ALL lanes (sorted by name)."""
    info = {}
    for p in sorted(Path(raw_dir).glob("*.fastq.gz")):
        m = ILLUMINA_FASTQ_RE.match(p.name) or SIMPLE_FASTQ_RE.match(p.name)
        if not m:
            continue
        entry = info.setdefault(m["sample"], {"R1": [], "R2": [], "type": "tumor", "sex": "", "purity": ""})
        entry[m["read"]].append(str(p))
    for sample in config.get("pon_samples") or []:
        if sample in info:
            info[sample]["type"] = "normal"
    return info

def load_samplesheet(path):
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#").fillna("")
    missing = {"sample", "fastq_1", "fastq_2"} - set(df.columns)
    if missing:
        raise ValueError(f"Sample sheet {path} is missing required column(s): {missing}")
    if df["sample"].duplicated().any():
        raise ValueError(f"Duplicated sample names in {path}: {df.loc[df['sample'].duplicated(), 'sample'].tolist()}")
    info = {}
    for _, row in df.iterrows():
        sample = row["sample"].strip()
        stype = row.get("type", "").strip().lower() or "tumor"
        if stype not in {"tumor", "normal"}:
            raise ValueError(f"Sample '{sample}': type must be 'tumor' or 'normal', got '{stype}'")
        sex = row.get("sex", "").strip().lower()
        if sex not in VALID_SEX:
            raise ValueError(f"Sample '{sample}': sex must be male/female or empty, got '{sex}'")
        purity = row.get("purity", "").strip()
        if purity and not (0 < float(purity) <= 1):
            raise ValueError(f"Sample '{sample}': purity must be in (0, 1], got '{purity}'")
        info[sample] = {
            "R1": [_resolve_fastq(f) for f in row["fastq_1"].split(",") if f.strip()],
            "R2": [_resolve_fastq(f) for f in row["fastq_2"].split(",") if f.strip()],
            "type": stype, "sex": VALID_SEX[sex], "purity": purity,
        }
    return info

SAMPLE_INFO = load_samplesheet(config["samplesheet"]) if config.get("samplesheet") else discover_samples(RAW_DIR)

# Sanity checks: paired files, same number of lanes, files exist
for _s, _i in SAMPLE_INFO.items():
    if not _i["R1"] or len(_i["R1"]) != len(_i["R2"]):
        raise ValueError(f"Sample '{_s}': R1/R2 FASTQ lists are empty or have a different number of lanes")
    _missing = [f for f in _i["R1"] + _i["R2"] if not Path(f).exists()]
    if _missing:
        raise FileNotFoundError(f"Sample '{_s}': FASTQ file(s) not found: {_missing}")

SAMPLES = sorted(SAMPLE_INFO)
NORMALS = [s for s in SAMPLES if SAMPLE_INFO[s]["type"] == "normal"]
TUMORS = [s for s in SAMPLES if SAMPLE_INFO[s]["type"] == "tumor"]

if not SAMPLES:
    raise ValueError(f"No samples found (sample sheet: {config.get('samplesheet')}, raw_dir: {RAW_DIR})")

print(f"Tumour samples: {TUMORS}")
print(f"Normal samples (CNV panel of normals): {NORMALS}")

# Restrict wildcards to the known sample names (avoids ambiguous matches such as
# sample='S1.filtered' when several rules write '<sample>.*.vcf.gz' in one folder)
wildcard_constraints:
    sample = "|".join(re.escape(s) for s in SAMPLES),
    normal = "|".join(re.escape(s) for s in SAMPLES),
    read = "R1|R2",
    caller = "GATK|VarDict"

# Set paths
QC_DIR = f"{RES_DIR}/01_RawReads_QC/02_Fastp"
ALIGN_DIR = f"{RES_DIR}/02_Alignment/Alignments"
ALIGN_QC_DIR = f"{RES_DIR}/02_Alignment/QC"
SMALL_DIR = f"{RES_DIR}/03_VariantCalling/01_SmallVariants"
GATK_DIR = f"{SMALL_DIR}/GATK"
VAR_DIR = f"{SMALL_DIR}/VarDict"
LUMPY_DIR = f"{RES_DIR}/03_VariantCalling/02_StructuralVariants/Lumpy"
MANTA_DIR = f"{RES_DIR}/03_VariantCalling/02_StructuralVariants/Manta"
CNV_DIR = f"{RES_DIR}/03_VariantCalling/03_CNVs"
GATK_CNV_DIR = f"{CNV_DIR}/GATK"
CNVkit_DIR = f"{CNV_DIR}/CNVkit"

# =============================================================================
# PIPELINE CONFIGURATION & VALIDATION
# =============================================================================
# Set type of variant calls to run based on config
VALID_TYPES = {"small", "sv", "cnv"}

# Allow user to specify a single type as a string, or multiple types as a list
calls_input = config["calls"]
if isinstance(calls_input, str):
    calls_input = [calls_input]

# Validate that the input is either 'all' or a list of valid types
if "all" in calls_input:
    CALL_TYPES = ["small", "sv", "cnv"]
else:
    # Validate that all specified call types are valid
    invalid_calls = set(calls_input) - VALID_TYPES
    if invalid_calls:
        raise ValueError(
            f"Invalid call type(s) specified: {invalid_calls}. "
            f"Choose from 'small', 'sv', 'cnv', or 'all'."
        )
    # Assign the validated list of call types
    CALL_TYPES = calls_input

CALLER_CHOICES = {
    "small": {"gatk": ["gatk"], "vardict": ["vardict"], "both": ["gatk", "vardict"]},
    "sv": {"manta": ["manta"], "lumpy": ["lumpy"], "both": ["manta", "lumpy"]},
    "cnv": {"cnvkit": ["cnvkit"], "gatk_cnv": ["gatk_cnv"], "both": ["cnvkit", "gatk_cnv"]},
}

def resolve_callers(kind):
    choice = config["callers"][kind]
    if choice not in CALLER_CHOICES[kind]:
        raise ValueError(f"Invalid {kind} caller '{choice}'. Choose from {list(CALLER_CHOICES[kind])}.")
    return CALLER_CHOICES[kind][choice]

SNP_CALLERS = resolve_callers("small")
SV_CALLERS = resolve_callers("sv")
CNV_CALLERS = resolve_callers("cnv")

# =============================================================================
# CNV PANEL OF NORMALS
# =============================================================================
CNV_CFG = config.get("cnv") or {}
PON_MODE = CNV_CFG.get("pon_mode", "build")
if PON_MODE not in ("build", "prebuilt"):
    raise ValueError(f"cnv.pon_mode must be 'build' or 'prebuilt', got '{PON_MODE}'")

if PON_MODE == "build":
    PON_DIR = f"{CNV_DIR}/00_PanelOfNormals"
else:
    PON_DIR = str(CNV_CFG.get("pon_dir") or "").rstrip("/")
    if "cnv" in CALL_TYPES and not Path(PON_DIR).is_dir():
        raise ValueError(f"cnv.pon_mode is 'prebuilt' but cnv.pon_dir does not exist: '{PON_DIR}'")

PON_QC = False
if "cnv" in CALL_TYPES and PON_MODE == "build":
    if len(NORMALS) < 2:
        raise ValueError(
            f"Building the CNV panel of normals requires at least 2 normal samples, found {len(NORMALS)}. "
            f"Mark them with type=normal in the sample sheet (or list them in 'pon_samples'), "
            f"or use cnv.pon_mode: prebuilt."
        )
    _recommended = int(CNV_CFG.get("pon_recommended_normals", 10))
    if len(NORMALS) < _recommended:
        print(f"[WARNING] Only {len(NORMALS)} normal samples for the CNV panel of normals "
              f"(recommended: >= {_recommended}). Expect noisier CNV calls.")
    PON_QC = bool(CNV_CFG.get("pon_qc", True))
    if PON_QC and len(NORMALS) < 3:
        print("[WARNING] Leave-one-out PoN QC needs >= 3 normal samples: skipped.")
        PON_QC = False

# Prebuilt PoN bundles are checked against the current BED before use
PON_CHECK = [] if PON_MODE == "build" else [f"{CNV_DIR}/Logs/pon_bundle_check.ok"]

# Sex chromosome contig names in the reference (excluded from GATK CNV if requested)
with open(f"{REFS['genome']}.fai") as _fai:
    SEX_CONTIGS = [l.split("\t")[0] for l in _fai if re.match(r"^(chr)?[XY]\t", l, re.IGNORECASE)]

# =============================================================================
# SHARED HELPERS
# =============================================================================
def civic_inputs(sv_dir):
    """Input function factory for the CIViC annotation of SV callers."""
    def _inputs(wildcards):
        inputs = {
            "tsv": f"{sv_dir}/02_Annotated/{wildcards.sample}_ann_gene.tsv",
            "civic": str(REFS['civic'])
        }
        gene_panel = config.get("fusion_gene_panel")
        if gene_panel:
            inputs["genes_list"] = gene_panel
        return inputs
    return _inputs

# =============================================================================
# DYNAMIC INPUT BUILDING FOR TARGET RULE
# =============================================================================
# Always include core QC and alignment reports
final_outputs = [
    f"{RES_DIR}/01_RawReads_QC/03_MultiQC/MultiQC_report.html",
    f"{RES_DIR}/02_Alignment/QC/multiqc_report.html"
]

# Append small variant workflow outputs if requested
if "small" in CALL_TYPES:
    if "gatk" in SNP_CALLERS:
        final_outputs.extend(expand(f"{GATK_DIR}/02_Annotated/{{sample}}_report.tsv", sample=TUMORS))
    if "vardict" in SNP_CALLERS:
        final_outputs.extend(expand(f"{VAR_DIR}/02_Annotated/{{sample}}_report.tsv", sample=TUMORS))

# Append structural variant workflow outputs if requested
if "sv" in CALL_TYPES:
    if "manta" in SV_CALLERS:
        final_outputs.extend(expand(f"{MANTA_DIR}/02_Annotated/{{sample}}_ann_civic.xlsx", sample=TUMORS))
    if "lumpy" in SV_CALLERS:
        final_outputs.extend(expand(f"{LUMPY_DIR}/02_Annotated/{{sample}}_ann_civic.xlsx", sample=TUMORS))

# Append copy number variant workflow outputs if requested
if "cnv" in CALL_TYPES:
    final_outputs.extend(expand(f"{CNV_DIR}/Report/{{sample}}_CNV_report.xlsx", sample=TUMORS))
    if "cnvkit" in CNV_CALLERS:
        final_outputs.extend(
            expand([
                f"{CNVkit_DIR}/04_Call/{{sample}}_scatter.pdf",
                f"{CNVkit_DIR}/04_Call/{{sample}}_diagram.pdf",
            ], sample=TUMORS)
        )
    if "gatk_cnv" in CNV_CALLERS:
        final_outputs.extend(
            expand([
                f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoised.png",
                f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.modeled.png",
            ], sample=TUMORS)
        )
    if PON_MODE == "build":
        final_outputs.append(f"{PON_DIR}/pon_manifest.tsv")
        if PON_QC:
            final_outputs.append(f"{PON_DIR}/pon_qc.tsv")

# =============================================================================
# TARGET RULE - Defines final expected outputs for the entire pipeline
# =============================================================================
rule all:
    input:
        final_outputs

# Build only the CNV panel of normals (e.g. from a sample sheet of healthy donors):
#   snakemake -s Snakefile.smk --configfile config.yml build_pon ...
rule build_pon:
    input:
        [f"{PON_DIR}/pon_manifest.tsv"] + ([f"{PON_DIR}/pon_qc.tsv"] if PON_QC else [])

# =============================================================================
# MODULE INCLUSIONS
# =============================================================================
include: SMK['qc']
include: SMK['ref']
if config["umi"]["strat"] == "none":
    include: SMK['align']
else:
    include: SMK['align_umi']
include: SMK['align_qc']
include: SMK['vc_common']
include: SMK['vc_gatk']
include: SMK['vc_vardict']
include: SMK['sv_manta']
include: SMK['sv_lumpy']
include: SMK['cnv_cnvkit']
include: SMK['cnv_gatk']
include: SMK['cnv_common']
