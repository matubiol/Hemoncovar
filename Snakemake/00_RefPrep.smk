"""
PIPELINE: Reference preparation
DESCRIPTION:
Writes the working copy of the intervals bed used by every downstream module
(Pre_processing/00_Harmonize_BED.py): chromosome naming harmonized with the
reference fasta, intervals on contigs absent from the reference removed,
header lines removed and a name (4th) column guaranteed.
"""

# =============================================================================
# PATHS
# =============================================================================
REF_DIR_RESULTS = f"{RES_DIR}/00_Reference"
INTERVALS = f"{REF_DIR_RESULTS}/{Path(RAW_BED).stem}_harmonized.bed"

# =============================================================================
# PIPELINE RULES
# =============================================================================
rule harmonize_bed_naming:
    input:
        bed = RAW_BED,
        fai = GENOME_FAI
    output:
        bed = INTERVALS
    params:
        needs_harmonization = NEEDS_BED_HARMONIZATION,
        ref_has_chr = fasta_has_chr_prefix(REFS['genome'])
    log:
        f"{REF_DIR_RESULTS}/Logs/harmonize_bed_naming.log"
    conda:
        ENVS['python']
    script:
        PRE['harmonize_bed']
