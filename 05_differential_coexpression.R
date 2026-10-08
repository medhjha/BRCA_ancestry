# Stage 5c: differential co-expression between EUR and AFR.
# First clean the gene universe to remove a low-expression artifact, then run
# DGCA on a pre-specified set of within-module and immune pairs.

library(DGCA)
library(data.table)
library(arrow)
library(yaml)
library(SummarizedExperiment)
library(matrixStats)

cfg <- read_yaml("config.yaml")
interim <- file.path(cfg$data, "interim")
proc <- file.path(cfg$data, "processed")
results <- cfg$results
set.seed(cfg$seeds$global)

coh <- as.data.table(read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
eur <- coh[patient12 %in% common & cohort == "EUR_primary", patient12]
afr <- coh[patient12 %in% common & cohort == "AFR_primary", patient12]

# clean the gene universe. the most variable genes include low-expression
# non-coding RNAs that give spurious near-perfect correlations resting on a few
# samples. keep protein-coding genes detected in at least half of samples in both
# groups, with a median count of at least 10, then re-rank by variance.
is_tumor <- function(bc) substr(bc, 14, 15) %in% sprintf("%02d", 1:9)
rna <- readRDS(file.path(interim, "tcga_rna_se.rds"))
rd <- as.data.table(as.data.frame(rowData(rna)))
rna_t <- rna[, is_tumor(colnames(rna))]
rna_t <- rna_t[, order(substr(colnames(rna_t), 14, 16))]
rna_t <- rna_t[, !duplicated(substr(colnames(rna_t), 1, 12))]
colnames(rna_t) <- substr(colnames(rna_t), 1, 12)
counts <- assay(rna_t, "unstranded")

protein_coding <- rd[gene_type == "protein_coding", gene_id]
g <- rownames(counts)[rownames(counts) %in% protein_coding]
det_eur <- rowMeans(counts[g, colnames(counts) %in% eur] > 0)
det_afr <- rowMeans(counts[g, colnames(counts) %in% afr] > 0)
g <- g[det_eur >= 0.5 & det_afr >= 0.5]
med_eur <- rowMedians(counts[g, colnames(counts) %in% eur])
med_afr <- rowMedians(counts[g, colnames(counts) %in% afr])
g <- g[med_eur >= 10 & med_afr >= 10]

vsd <- readRDS(file.path(proc, "rna_vst.rds"))
vsd_clean <- vsd[rownames(vsd) %in% g, ]
v <- rowVars(vsd_clean[, colnames(vsd_clean) %in% eur])
clean_genes <- rownames(vsd_clean)[order(v, decreasing = TRUE)[1:min(5000, nrow(vsd_clean))]]
saveRDS(clean_genes, file.path(proc, "dgca_genes_clean.rds"))

# pre-specify the pairs to test: all within-module pairs plus all pairs among a
# fixed immune gene panel. this is defined without looking at the result.
memb <- readRDS(file.path(proc, "module_membership.rds"))
within_pairs <- rbindlist(lapply(names(memb), function(m) {
  if (m == "0") return(NULL)
  gg <- intersect(memb[[m]], clean_genes)
  if (length(gg) < 2) return(NULL)
  cb <- t(combn(sort(gg), 2))
  data.table(Gene1 = cb[, 1], Gene2 = cb[, 2], set = paste0("module_", m))
}))
# add the pre-specified immune gene panel: all pairs among these genes that are
# in the cleaned universe. This is the panel used in the paper (immune-checkpoint,
# cytotoxic, B-cell, T-cell and myeloid markers).
immune_panel <- c(
  "CD19","MS4A1","CD79A","CD79B","CLEC17A","TCL1A","VPREB3","CD27",
  "CD3D","CD3E","CD3G","CD2","CD8A","CD8B","GZMA","GZMB","GZMK","PRF1","NKG7",
  "IFNG","CXCL9","CXCL10","CXCL11","CXCL13","CCL5","IL2RA","FOXP3","CTLA4",
  "PDCD1","CD274","LAG3","HAVCR2","TIGIT","ICOS","TNFRSF9",
  "CD68","CD163","CSF1R","ITGAX","LYZ","HLA-A","HLA-B","HLA-DRA","HLA-DRB1",
  "B2M","TAP1","TAP2")
immune_in <- intersect(immune_panel, clean_genes)
if (length(immune_in) >= 2) {
  ip <- t(combn(sort(immune_in), 2))
  immune_pairs <- data.table(Gene1 = ip[, 1], Gene2 = ip[, 2], set = "immune")
  pairset <- rbind(within_pairs, immune_pairs)
  pairset <- pairset[!duplicated(paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2)))]
} else {
  pairset <- within_pairs
}
pairset[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]
pair_genes <- unique(c(pairset$Gene1, pairset$Gene2))
saveRDS(pairset, file.path(proc, "dgca_pairset.rds"))
saveRDS(pair_genes, file.path(proc, "dgca_pairset_genes.rds"))

# run DGCA per module to keep memory bounded, then combine and do one global
# BH correction across the pre-specified pairs. this uses the parametric p-value;
# the permutation FDR is computed in 05b.
design_for <- function(samples) {
  grp <- ifelse(samples %in% eur, "EUR", "AFR")
  d <- model.matrix(~0 + factor(grp, levels = c("EUR", "AFR")))
  colnames(d) <- c("EUR", "AFR")
  d
}
run_module <- function(genes) {
  genes <- intersect(genes, clean_genes)
  if (length(genes) < 3) return(NULL)
  expr <- vsd[genes, colnames(vsd) %in% c(eur, afr)]
  design <- design_for(colnames(expr))
  dc <- ddcorAll(inputMat = expr, design = design, compare = c("EUR", "AFR"),
                 adjust = "none", nPerms = 0, corrType = "pearson", nPairs = "all")
  setDT(dc)
  dc
}
res <- rbindlist(lapply(setdiff(names(memb), "0"), function(m) run_module(memb[[m]])),
                 fill = TRUE)
res[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]
res <- res[!duplicated(key)][key %in% pairset$key]
res[, fdr := p.adjust(pValDiff, method = "BH")]

sym <- setNames(rowData(rna)$gene_name, rownames(rna))
res[, G1 := sym[Gene1]][, G2 := sym[Gene2]]
cat("pairs tested:", nrow(res), "\n")
cat("differential at FDR < 0.05:", sum(res$fdr < 0.05, na.rm = TRUE), "\n")

saveRDS(res, file.path(results, "dgca_final_EUR_vs_AFR.rds"))
fwrite(res[fdr < 0.05][order(fdr)], file.path(results, "dgca_significant_edges.csv"))
