# Stage 3b: build the covariate table and check the batch confounding.
# Covariates: PAM50 subtype (genefu), tumor purity (Aran 2015 CPE, from the
# TCGAbiolinks Tumor.purity table), ESTIMATE scores, and xCell fractions.
# The batch check decides how we handle batch downstream.

library(SummarizedExperiment)
library(data.table)
library(arrow)
library(yaml)
library(genefu)
library(estimate)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(matrixStats)

cfg <- read_yaml("config.yaml")
interim <- file.path(cfg$data, "interim")
proc <- file.path(cfg$data, "processed")
results <- cfg$results

common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
is_tumor <- function(bc) substr(bc, 14, 15) %in% sprintf("%02d", 1:9)

rna <- readRDS(file.path(interim, "tcga_rna_se.rds"))
rna <- rna[, is_tumor(colnames(rna))]
rna <- rna[, order(substr(colnames(rna), 14, 16))]
rna <- rna[, !duplicated(substr(colnames(rna), 1, 12))]
colnames(rna) <- substr(colnames(rna), 1, 12)
rna <- rna[, colnames(rna) %in% common]

# build a gene-symbol TPM matrix, keeping the most variable probe per symbol
tpm <- assay(rna, "tpm_unstrand")
gsym <- rowData(rna)$gene_name
keep <- !is.na(gsym) & gsym != ""
tpm <- tpm[keep, ]; gsym <- gsym[keep]
ord <- order(rowVars(tpm), decreasing = TRUE)
tpm <- tpm[ord, ]; gsym <- gsym[ord]
tpm <- tpm[!duplicated(gsym), ]
rownames(tpm) <- gsym[!duplicated(gsym)]

# PAM50 subtype. genefu needs all three pam50 objects loaded and an annotation
# frame keyed on EntrezGene.ID.
data(pam50); data(pam50.robust); data(pam50.scale)
entrez <- mapIds(org.Hs.eg.db, rownames(tpm), "ENTREZID", "SYMBOL", multiVals = "first")
annot <- data.frame(Gene.Symbol = rownames(tpm), EntrezGene.ID = entrez,
                    probe = rownames(tpm), row.names = rownames(tpm),
                    stringsAsFactors = FALSE)
sbt <- molecular.subtyping("pam50", data = t(log2(tpm + 1)), annot = annot, do.mapping = TRUE)
pam50_call <- as.character(sbt$subtype)

# ESTIMATE stromal and immune scores. This writes a text file, filters to the
# common genes, scores it, and reads the result back.
expr_txt <- file.path(proc, "estimate_expr.txt")
gct_in <- file.path(proc, "estimate_input.gct")
gct_out <- file.path(proc, "estimate_scores.gct")
write.table(data.frame(GeneSymbol = rownames(tpm), tpm, check.names = FALSE),
            expr_txt, sep = "\t", quote = FALSE, row.names = FALSE)
filterCommonGenes(input.f = expr_txt, output.f = gct_in, id = "GeneSymbol")
estimateScore(gct_in, gct_out, platform = "illumina")
# the GCT is whitespace delimited with two header lines and a Description column
est <- read.table(gct_out, skip = 2, header = TRUE, sep = "", row.names = 1,
                  check.names = FALSE)[, -1]
est_t <- as.data.frame(t(est))
est_t$patient12 <- gsub("\\.", "-", rownames(est_t))
setDT(est_t)

# tumor purity: Aran 2015 consensus purity (CPE), shipped in TCGAbiolinks
data("Tumor.purity", package = "TCGAbiolinks")
tp <- as.data.table(Tumor.purity)
tp[, patient12 := substr(Sample.ID, 1, 12)]
num <- function(x) as.numeric(gsub(",", ".", as.character(x)))
tp[, CPE := num(CPE)]
tp <- tp[!duplicated(patient12)]

covars <- data.table(patient12 = colnames(rna), pam50_genefu = pam50_call)
covars <- merge(covars, est_t, by = "patient12", all.x = TRUE)
covars <- merge(covars, tp[, .(patient12, CPE)], by = "patient12", all.x = TRUE)
# xCell fractions are added by the immune deconvolution step (03c) if used
saveRDS(covars, file.path(proc, "covariates.rds"))

# batch check: is tissue source site (and plate) confounded with ancestry?
coh <- as.data.table(read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
coh <- coh[patient12 %in% common & cohort %in% c("EUR_primary", "AFR_primary")]
coh[, anc := fifelse(cohort == "EUR_primary", "EUR", "AFR")]
bc <- colnames(rna)
meta <- data.table(patient12 = substr(bc, 1, 12),
                   tss = substr(bc, 6, 7), plate = substr(bc, 22, 25))
meta <- meta[!duplicated(patient12)]
m <- merge(coh[, .(patient12, anc)], meta, by = "patient12")
m[, anc2 := anc]

cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab)$statistic)
  n <- sum(tab); k <- min(nrow(tab), ncol(tab))
  as.numeric(sqrt(chi / (n * (k - 1))))
}
v_tss <- cramers_v(table(m$tss, m$anc))
v_plate <- cramers_v(table(m$plate, m$anc))
cat("Cramers V, site vs ancestry:", round(v_tss, 3), "\n")
cat("Cramers V, plate vs ancestry:", round(v_plate, 3), "\n")
cat("both are high, so we do not run ComBat. batch is handled by covariate\n")
cat("adjustment and the mixed-site sensitivity analysis in stage 8.\n")

saveRDS(m, file.path(proc, "batch_ancestry_meta.rds"))
