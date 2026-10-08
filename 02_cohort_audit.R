# Stage 2: find the patients that have all three layers, join to ancestry,
# and freeze the cohort. We use a 3-layer design (RNA, methylation, copy number)
# because requiring all four layers dropped the AFR sample size too far.

library(SummarizedExperiment)
library(data.table)
library(arrow)
library(yaml)

cfg <- read_yaml("config.yaml")
interim <- file.path(cfg$data, "interim")
proc <- file.path(cfg$data, "processed")
results <- cfg$results
dir.create(proc, showWarnings = FALSE, recursive = TRUE)

# helper: get the 12 character patient barcode from full sample barcodes
patient_id <- function(x) unique(substr(x, 1, 12))
is_tumor <- function(bc) substr(bc, 14, 15) %in% sprintf("%02d", 1:9)

rna <- readRDS(file.path(interim, "tcga_rna_se.rds"))
meth <- readRDS(file.path(interim, "tcga_meth_se.rds"))
cna <- readRDS(file.path(interim, "tcga_cna_se.rds"))

# copy number columns are paired "tumor;normal" barcodes, so split and take
# the tumor half
cna_barcodes <- unlist(strsplit(colnames(cna), ";"))
cna_tumor <- cna_barcodes[is_tumor(cna_barcodes)]

layers <- list(
  RNA = patient_id(colnames(rna)),
  METH = patient_id(colnames(meth)),
  CNA = patient_id(cna_tumor)
)
paired <- Reduce(intersect, layers)
cat("patients with all 3 layers:", length(paired), "\n")

# ancestry (parsed by python into a parquet keyed on patient barcode)
anc <- as.data.table(read_parquet(file.path(interim, "tcga_ancestry.parquet")))
anc <- anc[grepl("BRCA", toupper(tumor_type))]
anc[, patient12 := substr(patient_barcode, 1, 12)]

# assign each paired patient a role. strict EUR and strict AFR are the primary
# comparison; afr_admix is held out for the gradient analysis; everything else
# is excluded (too few, or admixed in a way we do not model here).
d <- anc[patient12 %in% paired]
d[, cohort := fifelse(consensus_ancestry == "eur", "EUR_primary",
              fifelse(consensus_ancestry == "afr", "AFR_primary",
              fifelse(consensus_ancestry == "afr_admix", "AFR_gradient", "exclude")))]

cat("\ncohort roles:\n")
print(d[, .N, by = cohort][order(-N)])

write_parquet(d, file.path(interim, "brca_cohort_3layer.parquet"))
fwrite(d, file.path(results, "brca_cohort_3layer.csv"))
cat("\nsaved cohort assignment\n")
