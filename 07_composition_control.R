# Stage 7: composition control and directional interpretation.
# DGCA cannot take covariates, so we residualize expression on purity and an
# independent adipocyte score, re-run the test, and keep the edges that survive
# (the clean-core set). Hub direction and enrichment are in 12_directional_enrichment.R.

library(DGCA)
library(xCell)
library(data.table)
library(arrow)
library(yaml)
library(SummarizedExperiment)

mem.maxVSize(48 * 1024)
cfg <- read_yaml("config.yaml")
interim <- file.path(cfg$data, "interim")
proc <- file.path(cfg$data, "processed")
results <- cfg$results
set.seed(cfg$seeds$global)

coh <- as.data.table(read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
eur <- coh[patient12 %in% common & cohort == "EUR_primary", patient12]
afr <- coh[patient12 %in% common & cohort == "AFR_primary", patient12]
samp <- c(eur, afr)

# independent adipocyte score from full xCell (not the 12-gene panel we test)
is_tumor <- function(bc) substr(bc, 14, 15) %in% sprintf("%02d", 1:9)
rna <- readRDS(file.path(interim, "tcga_rna_se.rds"))
rna_t <- rna[, is_tumor(colnames(rna))]
rna_t <- rna_t[, order(substr(colnames(rna_t), 14, 16))]
rna_t <- rna_t[, !duplicated(substr(colnames(rna_t), 1, 12))]
colnames(rna_t) <- substr(colnames(rna_t), 1, 12)
sym <- rowData(rna_t)$gene_name
tpm <- assay(rna_t, "tpm_unstrand")
keep <- !is.na(sym) & sym != ""
tpm <- tpm[keep, ]; rownames(tpm) <- make.unique(sym[keep])
tpm <- tpm[, colnames(tpm) %in% samp]
xc <- xCellAnalysis(tpm)
adipo <- xc["Adipocytes", ]

cov <- as.data.table(readRDS(file.path(proc, "covariates.rds")))
purity <- setNames(cov$CPE, cov$patient12)

vsd <- readRDS(file.path(proc, "rna_vst.rds"))
genes <- readRDS(file.path(proc, "dgca_genes_clean.rds"))
E <- vsd[genes, samp]
pur_v <- purity[colnames(E)]; adi_v <- adipo[colnames(E)]
ok <- !is.na(pur_v) & !is.na(adi_v)
E <- E[, ok]; grp <- ifelse(colnames(E) %in% eur, "EUR", "AFR")
pur_v <- pur_v[ok]; adi_v <- adi_v[ok]

# residualize each gene on purity and the adipocyte score
resid_on <- function(mat, covs) {
  X <- model.matrix(~., data = as.data.frame(covs))
  t(residuals(lm(t(mat) ~ X - 1)))
}
E_adj <- resid_on(E, data.frame(purity = pur_v, adipo = adi_v))

memb <- readRDS(file.path(proc, "module_membership.rds"))
pairset <- as.data.table(readRDS(file.path(proc, "dgca_pairset.rds")))
pairset[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]
design <- model.matrix(~0 + factor(grp, levels = c("EUR", "AFR")))
colnames(design) <- c("EUR", "AFR")

run_adj <- function(mat) {
  res <- rbindlist(lapply(setdiff(names(memb), "0"), function(m) {
    gg <- intersect(memb[[m]], rownames(mat))
    if (length(gg) < 3) return(NULL)
    dc <- ddcorAll(inputMat = mat[gg, ], design = design, compare = c("EUR", "AFR"),
                   adjust = "none", nPerms = 0, corrType = "pearson", nPairs = "all")
    setDT(dc)
    dc
  }), fill = TRUE)
  res[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]
  res <- res[!duplicated(key)][key %in% pairset$key]
  res[, fdr_pa := p.adjust(pValDiff, method = "BH")]
  res
}
adj <- run_adj(E_adj)

full <- as.data.table(readRDS(file.path(results, "dgca_final_EUR_vs_AFR.rds")))
adipo_panel <- c("PLIN1", "PLIN4", "RBP4", "GPD1", "AQP7", "PPARG", "ADIPOQ",
                 "LEP", "FABP4", "CD36", "LPL", "CFD")
cmp <- merge(full[, .(key, G1, G2, fdr_orig = fdr)],
             adj[, .(key, fdr_pa)], by = "key", all.x = TRUE)
cmp[, is_adipo_edge := (G1 %in% adipo_panel | G2 %in% adipo_panel)]
saveRDS(cmp, file.path(results, "dgca_adjusted_comparison.rds"))
cat("edges surviving composition adjustment (clean-core):",
    cmp[fdr_orig < 0.05 & fdr_pa < 0.05, .N], "\n")
