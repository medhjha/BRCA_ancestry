# Stage 5c permutation FDR.
# The parametric FDR in 05 is fast but a permutation FDR is the stricter test we
# report. DGCA's splitSet builds the permutation null once across the full matrix,
# so we split the genes into small batches, run each batch against the full
# matrix, and combine. splitSet returns pairs linking the split genes to the
# rest, so we use two offset passes plus a targeted fill to cover the within-batch
# pairs too. This is memory-safe on a machine with limited RAM.

library(DGCA)
library(data.table)
library(arrow)
library(yaml)
library(SummarizedExperiment)

mem.maxVSize(20 * 1024)
cfg <- read_yaml("config.yaml")
interim <- file.path(cfg$data, "interim")
proc <- file.path(cfg$data, "processed")
results <- cfg$results
set.seed(cfg$seeds$global)

coh <- as.data.table(read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
eur <- coh[patient12 %in% common & cohort == "EUR_primary", patient12]
afr <- coh[patient12 %in% common & cohort == "AFR_primary", patient12]
vsd <- readRDS(file.path(proc, "rna_vst.rds"))
pair_genes <- readRDS(file.path(proc, "dgca_pairset_genes.rds"))
pairset <- as.data.table(readRDS(file.path(proc, "dgca_pairset.rds")))
pairset[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]

expr <- vsd[pair_genes, colnames(vsd) %in% c(eur, afr)]
grp <- ifelse(colnames(expr) %in% eur, "EUR", "AFR")
design <- model.matrix(~0 + factor(grp, levels = c("EUR", "AFR")))
colnames(design) <- c("EUR", "AFR")

# run one pass over a given gene ordering in chunks, collecting the empirical
# p-values for the pre-specified pairs
run_pass <- function(gene_order, chunk_size = 100) {
  chunks <- split(gene_order, ceiling(seq_along(gene_order) / chunk_size))
  acc <- list()
  for (i in seq_along(chunks)) {
    s <- as.data.table(ddcorAll(inputMat = expr, design = design,
                                compare = c("EUR", "AFR"), adjust = "perm",
                                nPerms = 1000, corrType = "pearson",
                                nPairs = "all", splitSet = chunks[[i]]))
    s[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]
    acc[[i]] <- s[key %in% pairset$key, .(key, empPVals)]
    rm(s); gc(verbose = FALSE)
  }
  rbindlist(acc)
}

# two offset passes so within-chunk pairs of one pass become cross-chunk in the
# other, then a targeted fill for anything still missing
p1 <- run_pass(pair_genes)
shifted <- c(pair_genes[51:length(pair_genes)], pair_genes[1:50])
p2 <- run_pass(shifted)
res <- rbind(p1, p2)[!duplicated(key)]

missing <- pairset[!key %in% res$key]
if (nrow(missing) > 0) {
  fill_genes <- unique(missing$Gene1)
  p3 <- run_pass(fill_genes, chunk_size = 60)
  res <- rbind(res, p3[key %in% missing$key])[!duplicated(key)]
}
cat("coverage:", round(100 * nrow(res) / nrow(pairset), 1), "percent of pairs\n")

res[, fdr_perm := p.adjust(empPVals, method = "BH")]

# attach the composition classification (from stage 7) so we can pick out the
# clean-core set
cmp <- as.data.table(readRDS(file.path(results, "dgca_adjusted_comparison.rds")))
res <- merge(res, cmp[, .(key, fdr_pa, is_adipo_edge)], by = "key", all.x = TRUE)
res <- merge(res, pairset[, .(key, Gene1, Gene2)], by = "key", all.x = TRUE)
sym <- setNames(rowData(readRDS(file.path(interim, "tcga_rna_se.rds")))$gene_name,
                rownames(readRDS(file.path(interim, "tcga_rna_se.rds"))))
res[, G1 := sym[Gene1]][, G2 := sym[Gene2]]

clean_core <- res[fdr_perm < 0.05 & fdr_pa < 0.05]
cat("permutation FDR < 0.05:", sum(res$fdr_perm < 0.05, na.rm = TRUE), "\n")
cat("clean-core (composition robust):", nrow(clean_core), "\n")

saveRDS(res, file.path(results, "dgca_permFDR_final.rds"))
fwrite(res[fdr_perm < 0.10][order(fdr_perm)], file.path(results, "dgca_permFDR_significant.csv"))
