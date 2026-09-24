#!/usr/bin/env Rscript
# ============================================================================ #
# Clinical annotation of SVs / fusions with CIViC                               #
#   - Fusions ("GeneA::GeneB") are matched in both orientations.                 #
#   - Non-fusion rows are matched as "<Gene> DEL|DUP|INS".                        #
#   - Optional filter by a list of genes of interest (e.g. the 26 HS2 fusion     #
#     genes).                                                                    #
# Packages must be provided by the execution environment (conda/module): they  #
# are NOT installed at runtime.                                                 #
# ============================================================================ #
suppressPackageStartupMessages({
  for (p in c("dplyr", "tidyr", "optparse", "writexl")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("R package '%s' is not available in this environment.", p))
    }
    library(p, character.only = TRUE)
  }
})

# ============================================================================ #
#                         Setup command line arguments                         #
# ============================================================================ #
if (exists("snakemake")) {
  opt <- list(
    tsv = snakemake@input[["tsv"]],
    civic = snakemake@input[["civic"]],
    genes_list = snakemake@input[["genes_list"]],
    xlsx = snakemake@output[["xlsx"]]
  )
  if (length(snakemake@log) > 0) {
    log_con <- file(snakemake@log[[1]], open = "wt")
    sink(log_con, type = "message")
  }
} else {
  option_list <- list(
    make_option(c("-i", "--tsv"), type = "character"),
    make_option(c("-c", "--civic"), type = "character"),
    make_option(c("-g", "--genes_list"), type = "character", default = NULL,
                help = "Optional: file with one gene symbol per line to filter fusions/SVs by"),
    make_option(c("-o", "--xlsx"), type = "character")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
}

vars <- read.delim(opt$tsv, stringsAsFactors = FALSE, colClasses = "character")
message(sprintf("SV rows read: %d", nrow(vars)))

# ============================================================================ #
#                               Prepare CIViC                                  #
# ============================================================================ #
civic <- read.delim(opt$civic, stringsAsFactors = FALSE) %>%
  separate(
    col = molecular_profile,
    sep = " ",
    into = c("Var", "Var_type"),
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(
    Var_type = coalesce(Var_type, ""),
    Var_type_short = case_when(
      grepl("delins", Var_type, ignore.case = TRUE) ~ "COMPLEX",
      Var_type == "Deleterious" ~ "Harmful",
      grepl("del", Var_type, ignore.case = TRUE) ~ "DEL",
      grepl("ins", Var_type, ignore.case = TRUE) ~ "INS",
      grepl("Amplification", Var_type, ignore.case = TRUE) | grepl("dup", Var_type, ignore.case = TRUE) ~ "DUP",
      TRUE ~ Var_type
    ),
    Gene = case_when(
      Var_type == "Fusion" ~ Var,
      Var_type_short %in% c("DUP", "DEL", "INS") ~ paste(Var, Var_type_short),
      TRUE ~ Var
    )
  ) %>%
  select(Gene, Var_type, disease)

civic_fusions_reversed <- civic %>%
  filter(Var_type == "Fusion", grepl("::", Gene)) %>%
  separate(Gene, into = c("g1", "g2"), sep = "::", remove = FALSE) %>%
  mutate(Gene = paste0(g2, "::", g1)) %>%
  select(Gene, Var_type, disease)

civic <- bind_rows(civic, civic_fusions_reversed) %>%
  distinct()

# ============================================================================ #
#             Optional filter by list of genes of interest                     #
# ============================================================================ #
has_gene_filter <- !is.null(opt$genes_list) && length(opt$genes_list) > 0 &&
  !is.na(opt$genes_list) && nzchar(opt$genes_list)

vars_ann <- vars %>%
  separate(col = "Gene", into = c("Gene1", "Gene2"), sep = "::", remove = FALSE, fill = "right")

if (has_gene_filter) {
  genes_fus <- trimws(readLines(opt$genes_list))
  genes_fus <- genes_fus[nzchar(genes_fus)]
  message(sprintf("Gene filter enabled: %d genes loaded from %s", length(genes_fus), opt$genes_list))
  vars_ann <- vars_ann %>%
    filter(Gene1 %in% genes_fus | Gene2 %in% genes_fus)
} else {
  message("No genes_list provided: skipping gene-based filtering, keeping all variants.")
}

# ============================================================================ #
#                             Join with CIViC                                  #
# ============================================================================ #
# Fusions keep their "A::B" key; single-gene DEL/DUP/INS rows get the SV type
vars_ann <- vars_ann %>%
  select(-Gene1, -Gene2) %>%
  mutate(
    Join_key = if_else(
      !grepl("::", Gene) & SVType %in% c("DUP", "DEL", "INS"),
      paste(Gene, SVType),
      Gene
    )
  ) %>%
  left_join(civic, by = c("Join_key" = "Gene"), relationship = "many-to-many") %>%
  mutate(disease = if_else(disease == "", NA_character_, disease))

group_cols <- setdiff(colnames(vars_ann), c("Join_key", "Var_type", "disease"))

vars_ann <- vars_ann %>%
  group_by(across(all_of(group_cols))) %>%
  summarise(
    CIViC_disease = {
      d <- unique(na.omit(disease))
      if (length(d) == 0) NA_character_ else paste(d, collapse = ", ")
    },
    .groups = "drop"
  ) %>%
  mutate(Support = as.numeric(Support)) %>%
  arrange(desc(Fusion_candidate == "True"), desc(Support), CIViC_disease)

message(sprintf("Rows written: %d", nrow(vars_ann)))
write_xlsx(vars_ann, opt$xlsx)
