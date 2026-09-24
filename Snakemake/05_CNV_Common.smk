"""
PIPELINE: CNV panel of normals bundle, PoN QC and per-sample CNV report
DESCRIPTION:
  - Build mode: writes pon_manifest.tsv (normals, sexes, pipeline/wrapper
    versions, md5 of the intervals BED and of every PoN file) next to the
    CNVkit/GATK PoN files, so that the whole 00_PanelOfNormals folder can be
    reused later with cnv.pon_mode: prebuilt.
  - Prebuilt mode: checks that the bundle was built with the same BED.
  - PoN QC (leave-one-out): noise, false CNV calls, declared vs inferred sex
    and possible mosaic loss of chrY for every normal
    (Post_processing/05_CNV_PoN_QC.py).
  - Per tumour sample: one Excel report combining both callers
    (Post_processing/06_CNV_Report.py): arm-level events, altered segments,
    report genes, allelic imbalance / copy-neutral LOH and parameters.
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
REPORT_GENES = f"{RES_DIR}/00_Reference/cnv_report_genes.tsv"

def pon_files():
    """PoN files expected in the bundle for the selected CNV callers."""
    files = []
    if "cnvkit" in CNV_CALLERS:
        files += [f"{PON_DIR}/cnvkit/{f}" for f in ("targets.bed", "antitargets.bed", "reference.cnn")]
    if "gatk_cnv" in CNV_CALLERS:
        files += [f"{PON_DIR}/gatk/{f}" for f in ("targets.preprocessed.interval_list", "panel_of_normals.hdf5")]
    return files

# =============================================================================
# REPORT GENES
# =============================================================================
rule cnv_report_genes:
    input:
        bed = INTERVALS,
        gtf = REFS['gtf'],
        genes = [CNV_CFG["report_genes"]] if CNV_CFG.get("report_genes") else []
    output:
        REPORT_GENES
    log:
        f"{RES_DIR}/00_Reference/Logs/cnv_report_genes.log"
    conda:
        ENVS['python']
    script:
        PRE['report_genes']

# =============================================================================
# PANEL OF NORMALS BUNDLE
# =============================================================================
if PON_MODE == "build":
    rule pon_manifest:
        input:
            bed = INTERVALS,
            pon = pon_files()
        output:
            f"{PON_DIR}/pon_manifest.tsv"
        params:
            normals = ",".join(NORMALS),
            sexes = ",".join(f"{n}:{SAMPLE_INFO[n]['sex'] or 'unknown'}" for n in NORMALS),
            callers = ",".join(CNV_CALLERS),
            umi = config["umi"]["strat"],
            version = PIPELINE_VERSION,
            wrapper = WRAPPER,
            genome = str(GENOME)
        shell:
            """
            {{
            printf "key\\tvalue\\n"
            printf "created\\t%s\\n" "$(date -Iseconds)"
            printf "pipeline_version\\t{params.version}\\n"
            printf "wrapper_version\\t{params.wrapper}\\n"
            printf "genome\\t{params.genome}\\n"
            printf "umi_strategy\\t{params.umi}\\n"
            printf "callers\\t{params.callers}\\n"
            printf "n_normals\\t%s\\n" "$(echo {params.normals} | tr ',' '\\n' | wc -l)"
            printf "normals\\t{params.normals}\\n"
            printf "normal_sexes\\t{params.sexes}\\n"
            printf "intervals_md5\\t%s\\n" "$(md5sum {input.bed} | cut -d' ' -f1)"
            for f in {input.pon}; do
                printf "md5:%s\\t%s\\n" "$(basename $(dirname $f))/$(basename $f)" "$(md5sum $f | cut -d' ' -f1)"
            done
            }} > {output}
            """
else:
    # The bundle must have been built with the same (harmonized) intervals BED
    rule check_pon_bundle:
        input:
            bed = INTERVALS,
            manifest = f"{PON_DIR}/pon_manifest.tsv",
            pon = pon_files()
        output:
            f"{CNV_DIR}/Logs/pon_bundle_check.ok"
        shell:
            """
            expected=$(awk -F'\\t' '$1=="intervals_md5"{{print $2}}' {input.manifest})
            actual=$(md5sum {input.bed} | cut -d' ' -f1)
            if [ "$expected" != "$actual" ]; then
                echo "ERROR: the prebuilt PoN in {PON_DIR} was built with a different intervals BED" >&2
                echo "       (manifest md5: $expected, current BED md5: $actual)." >&2
                exit 1
            fi
            cp {input.manifest} {output}
            """

# =============================================================================
# PoN QC (build mode, leave-one-out)
# =============================================================================
if PON_QC:
    rule pon_qc:
        input:
            cnvkit_cnr = expand(f"{CNVkit_DIR}/PoN_QC/{{normal}}.cnr", normal=NORMALS) if "cnvkit" in CNV_CALLERS else [],
            cnvkit_call = expand(f"{CNVkit_DIR}/PoN_QC/{{normal}}.call.cns", normal=NORMALS) if "cnvkit" in CNV_CALLERS else [],
            cnvkit_sex = [f"{CNVkit_DIR}/PoN_QC/cnvkit_sex.tsv"] if "cnvkit" in CNV_CALLERS else [],
            gatk_dcr = expand(f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.denoisedCR.tsv", normal=NORMALS) if "gatk_cnv" in CNV_CALLERS else [],
            gatk_called = expand(f"{GATK_CNV_DIR}/PoN_QC/{{normal}}.called.seg", normal=NORMALS) if "gatk_cnv" in CNV_CALLERS else []
        output:
            f"{PON_DIR}/pon_qc.tsv"
        params:
            normals = NORMALS,
            sexes = [SAMPLE_INFO[n]["sex"] for n in NORMALS],
            min_probes = CNV_CFG.get("min_probes", 10)
        log:
            f"{PON_DIR}/Logs/pon_qc.log"
        conda:
            ENVS['python']
        script:
            SRC['CNV_pon_qc']

# =============================================================================
# PER-SAMPLE CNV REPORT
# =============================================================================
def cnv_report_inputs(wildcards):
    inputs = {"cyto": str(REFS['cyto']), "genes": REPORT_GENES}
    if "cnvkit" in CNV_CALLERS:
        inputs["cnvkit"] = f"{CNVkit_DIR}/05_Annotated/{wildcards.sample}.ann.tsv"
    if "gatk_cnv" in CNV_CALLERS:
        inputs["gatk"] = f"{GATK_CNV_DIR}/04_Annotated/{wildcards.sample}.ann.tsv"
        inputs["gatk_model"] = f"{GATK_CNV_DIR}/02_Segmentation/{wildcards.sample}.modelFinal.seg"
    return inputs

rule cnv_report:
    input:
        unpack(cnv_report_inputs)
    output:
        f"{CNV_DIR}/Report/{{sample}}_CNV_report.xlsx"
    params:
        sample = "{sample}",
        sex = lambda wildcards: SAMPLE_INFO[wildcards.sample]["sex"],
        purity = lambda wildcards: SAMPLE_INFO[wildcards.sample]["purity"],
        min_probes = CNV_CFG.get("min_probes", 10),
        arm_min_fraction = CNV_CFG.get("arm_min_fraction", 0.5),
        loh_max_maf = CNV_CFG.get("loh_max_maf", 0.35),
        loh_min_hets = CNV_CFG.get("loh_min_hets", 10),
        neutral_bounds = CNV_CFG.get("neutral_bounds", [0.9, 1.1]),
        pon_mode = PON_MODE,
        pon_dir = PON_DIR,
        version = PIPELINE_VERSION
    log:
        f"{CNV_DIR}/Report/Logs/{{sample}}.log"
    conda:
        ENVS['python']
    script:
        SRC['CNV_report']
