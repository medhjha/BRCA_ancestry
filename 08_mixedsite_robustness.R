# Stage 8: mixed-site robustness.
# Site is confounded with ancestry, so we re-run the differential analysis using
# only patients from sites that hold both EUR and AFR patients, and check whether
# the effect sizes agree with the full analysis. If they do, the edges are
# ancestry-driven, not a site artifact.

library(DGCA)
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

common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
coh <- as.data.table(read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
coh <- coh[patient12 %in% common & cohort %in% c("EUR_primary", "AFR_primary")]
coh[, anc := fifelse(cohort == "EUR_primary", "EUR", "AFR")]

# recover site from the RNA barcodes and find mixed sites
rna <- readRDS(file.path(interim, "tcga_rna_se.rds"))
bc <- data.table(patient12 = substr(colnames(rna), 1, 12), tss = substr(colnames(rna), 6, 7))
bc <- bc[!duplicated(patient12)]
m <- merge(coh[, .(patient12, anc)], bc, by = "patient12")
site_comp <- m[, .(nEUR = sum(anc == "EUR"), nAFR = sum(anc == "AFR")), by = tss]
mixed_sites <- site_comp[nEUR > 0 & nAFR > 0, tss]
mixed <- m[tss %in% mixed_sites]
eur <- mixed[anc == "EUR", patient12]
afr <- mixed[anc == "AFR", patient12]
cat("mixed-site cohort: EUR", length(eur), "AFR", length(afr), "\n")

vsd <- readRDS(file.path(proc, "rna_vst.rds"))
memb <- readRDS(file.path(proc, "module_membership.rds"))
clean_genes <- readRDS(file.path(proc, "dgca_genes_clean.rds"))
pairset <- as.data.table(readRDS(file.path(proc, "dgca_pairset.rds")))
pairset[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]

design <- model.matrix(~0 + factor(ifelse(c(eur, afr) %in% eur, "EUR", "AFR"),
                                   levels = c("EUR", "AFR")))
colnames(design) <- c("EUR", "AFR")
run_module <- function(genes) {
  genes <- intersect(genes, clean_genes)
  if (length(genes) < 3) return(NULL)
  expr <- vsd[genes, c(eur, afr)]
  dc <- ddcorAll(inputMat = expr, design = design, compare = c("EUR", "AFR"),
                 adjust = "none", nPerms = 0, corrType = "pearson", nPairs = "all")
  setDT(dc)
  dc
}
res <- rbindlist(lapply(setdiff(names(memb), "0"), function(m) run_module(memb[[m]])),
                 fill = TRUE)
res[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]
res <- res[!duplicated(key)][key %in% pairset$key]

# compare effect sizes with the full analysis
full <- as.data.table(readRDS(file.path(results, "dgca_final_EUR_vs_AFR.rds")))
full[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]
cmp <- merge(full[, .(key, z_full = zScoreDiff, fdr_full = fdr)],
             res[, .(key, z_mixed = zScoreDiff)], by = "key")
sp <- cor(cmp$z_full, cmp$z_mixed, method = "spearman", use = "complete.obs")
sig <- cmp[fdr_full < 0.05]
dir_agree <- mean(sign(sig$z_full) == sign(sig$z_mixed), na.rm = TRUE)
cat("effect-size Spearman, full vs mixed-site:", round(sp, 3), "\n")
cat("direction agreement on significant edges:", round(100 * dir_agree, 1), "percent\n")

fwrite(cmp, file.path(results, "dgca_full_vs_mixedsite.csv"))
