# Stage 3: preprocess the three layers and build the covariate table.
# RNA gets variance-stabilized, methylation goes to M-values, copy number takes
# the tumor half of the ASCAT pairs. Covariates are PAM50 subtype, tumor purity,
# ESTIMATE scores, and xCell fractions.

library(SummarizedExperiment)
library(DESeq2)
library(data.table)
library(arrow)
library(yaml)
library(matrixStats)

cfg <- read_yaml("config.yaml")
interim <- file.path(cfg$data, "interim")
proc <- file.path(cfg$data, "processed")

coh <- as.data.table(read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
coh <- coh[cohort != "exclude"]
cohort_patients <- coh$patient12

is_tumor <- function(bc) substr(bc, 14, 15) %in% sprintf("%02d", 1:9)

# keep one primary tumor sample per patient, then rename columns to the 12
# character patient barcode
dedup <- function(se) {
  se <- se[, is_tumor(colnames(se))]
  se <- se[, order(substr(colnames(se), 14, 16))]
  se <- se[, !duplicated(substr(colnames(se), 1, 12))]
  colnames(se) <- substr(colnames(se), 1, 12)
  se
}

# RNA: variance stabilizing transform, fit on all tumors then subset to cohort
rna <- dedup(readRDS(file.path(interim, "tcga_rna_se.rds")))
counts <- assay(rna, "unstranded")
storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
counts[is.na(counts)] <- 0L
col_data <- DataFrame(dummy = factor(rep(1, ncol(counts))))
rownames(col_data) <- colnames(counts)
dds <- DESeqDataSetFromMatrix(counts, colData = col_data, design = ~1)
vsd <- assay(vst(dds, blind = TRUE))
vsd <- vsd[, colnames(vsd) %in% cohort_patients, drop = FALSE]
saveRDS(vsd, file.path(proc, "rna_vst.rds"))

# Methylation: keep beta for interpretation, convert to M-values for modeling.
# Drop probes with more than 20 percent missing across the cohort.
meth <- dedup(readRDS(file.path(interim, "tcga_meth_se.rds")))
beta <- assay(meth)
beta <- beta[, colnames(beta) %in% cohort_patients, drop = FALSE]
beta <- beta[rowMeans(is.na(beta)) <= 0.20, , drop = FALSE]
eps <- 1e-3
beta_clamped <- pmin(pmax(beta, eps), 1 - eps)
mval <- log2(beta_clamped / (1 - beta_clamped))
saveRDS(beta, file.path(proc, "meth_beta.rds"))
saveRDS(mval, file.path(proc, "meth_mval.rds"))

# Copy number: take the tumor half of each pair, winsorize the extreme
# amplifications at the 99.9th percentile
cna <- readRDS(file.path(interim, "tcga_cna_se.rds"))
cn <- assay(cna, "copy_number")
pick_tumor <- function(colnm) {
  parts <- strsplit(colnm, ";")[[1]]
  tum <- parts[is_tumor(parts)]
  if (length(tum)) tum[1] else NA_character_
}
tum <- vapply(colnames(cn), pick_tumor, character(1))
cn <- cn[, !is.na(tum), drop = FALSE]
tum <- tum[!is.na(tum)]
colnames(cn) <- substr(tum, 1, 12)
cn <- cn[, !duplicated(colnames(cn)), drop = FALSE]
cn <- cn[, colnames(cn) %in% cohort_patients, drop = FALSE]
cap <- quantile(cn, 0.999, na.rm = TRUE)
cn[cn > cap] <- cap
saveRDS(cn, file.path(proc, "cna_tumor.rds"))

# the aligned cohort is the set of patients present in all three processed matrices
common <- Reduce(intersect, list(colnames(vsd), colnames(beta), colnames(cn)))
saveRDS(common, file.path(proc, "cohort_patients_aligned.rds"))
cat("aligned cohort:", length(common), "patients\n")
cat("preprocessing done. covariates are built in 03b_covariates.R\n")
