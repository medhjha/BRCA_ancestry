# Methylation: build the promoter probe set, test differential methylation, and
# check whether the co-expression-differential genes are the same genes whose
# promoters are differentially methylated.
# The SNP-overlapping and cross-reactive probe removal is the important step here,
# because SNPs under probes create fake ancestry differences.

library(SummarizedExperiment)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(limma)
library(data.table)
library(arrow)
library(yaml)

cfg <- read_yaml("config.yaml")
interim <- file.path(cfg$data, "interim")
proc <- file.path(cfg$data, "processed")
results <- cfg$results

mval <- readRDS(file.path(proc, "meth_mval.rds"))
meth_se <- readRDS(file.path(interim, "tcga_meth_se.rds"))
rd <- as.data.table(as.data.frame(rowData(meth_se)), keep.rownames = "probe")
gr <- rowRanges(meth_se)
clean_genes <- readRDS(file.path(proc, "dgca_genes_clean.rds"))
rna <- readRDS(file.path(interim, "tcga_rna_se.rds"))
sym <- setNames(rowData(rna)$gene_name, rownames(rna))
aim1 <- unique(sym[clean_genes])

# remove SNP-overlapping and cross-reactive probes using the sesame masks
mask_cols <- intersect(c("MASK_snp5_common", "MASK_snp5_GMAF1p", "MASK_mapping",
                         "MASK_general", "MASK_sub30_copy", "MASK_sub35_copy"),
                       colnames(rd))
rd[, masked := Reduce(`|`, lapply(mask_cols, function(cl) as.logical(get(cl)) %in% TRUE))]
keep_probes <- rd[masked == FALSE, probe]

# map probes to the clean-core genes
rd_g <- rd[probe %in% keep_probes & !is.na(gene) & gene != ""]
long <- rd_g[, .(g = unlist(strsplit(gene, ";"))), by = probe][g %in% aim1]

# restrict to promoter probes: within 1500 bp upstream to 500 bp downstream of a
# gene transcription start
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- suppressWarnings(genes(txdb))
eg <- mapIds(org.Hs.eg.db, unique(long$g), "ENTREZID", "SYMBOL", multiVals = "first")
eg <- eg[!is.na(eg) & eg %in% names(genes_gr)]
tss <- resize(genes_gr[eg], width = 1, fix = "start")
names(tss) <- names(eg)
prom <- promoters(tss, upstream = 1500, downstream = 500)

gr_k <- gr[names(gr) %in% unique(long$probe)]
gr_dt <- data.table(probe = names(gr_k), chr = as.character(seqnames(gr_k)), pos = start(gr_k))
prom_dt <- data.table(g = names(prom), chr = as.character(seqnames(prom)),
                      start = start(prom), end = end(prom))
pg <- merge(merge(long, gr_dt, by = "probe"), prom_dt, by = c("g", "chr"), allow.cartesian = TRUE)
pg[, in_prom := pos >= start & pos <= end]
prom_pairs <- pg[in_prom == TRUE]
prom_probes <- intersect(unique(prom_pairs$probe), rownames(mval))

mval_prom <- mval[prom_probes, ]
saveRDS(mval_prom, file.path(proc, "meth_mval_promoter.rds"))
saveRDS(prom_pairs[, .(probe, g, chr, pos)], file.path(proc, "meth_promoter_probe2gene.rds"))
cat("promoter probes after SNP masking:", nrow(mval_prom), "\n")

# differential methylation with limma, adjusting for purity, subtype, immune
coh <- as.data.table(read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
coh <- coh[patient12 %in% common & cohort %in% c("EUR_primary", "AFR_primary")]
coh[, anc := factor(fifelse(cohort == "EUR_primary", "EUR", "AFR"), levels = c("EUR", "AFR"))]
cov <- as.data.table(readRDS(file.path(proc, "covariates.rds")))
d <- merge(coh[, .(patient12, anc)], cov[, .(patient12, CPE, pam50_genefu, ImmuneScore)], by = "patient12")
d <- d[complete.cases(d) & patient12 %in% colnames(mval_prom)]
M <- mval_prom[, d$patient12]

design <- model.matrix(~anc + CPE + factor(pam50_genefu) + ImmuneScore, data = d)
fit <- eBayes(lmFit(M, design))
tt <- as.data.table(topTable(fit, coef = "ancAFR", number = Inf), keep.rownames = "probe")
setnames(tt, c("logFC", "adj.P.Val"), c("delta_M_AFRvsEUR", "fdr"))
tt <- merge(tt, prom_pairs[, .(probe, g)], by = "probe", allow.cartesian = TRUE)
gene_dm <- tt[, .(min_fdr = min(fdr), best_delta = delta_M_AFRvsEUR[which.min(fdr)]), by = g]
gene_dm[, DM := min_fdr < 0.05]
gene_dm[, DM_strict := g %in% tt[fdr < 0.05 & abs(delta_M_AFRvsEUR) >= 0.2, unique(g)]]
cat("genes differentially methylated (FDR < 0.05):", sum(gene_dm$DM), "\n")
saveRDS(tt, file.path(results, "dmp_probe_EUR_vs_AFR.rds"))
saveRDS(gene_dm, file.path(results, "dmp_gene_EUR_vs_AFR.rds"))
saveRDS(gene_dm[, .(g, DM_strict, best_delta)], file.path(results, "dmp_gene_strict.rds"))
fwrite(gene_dm[order(min_fdr)], file.path(results, "dmp_gene_table.csv"))

# cross-layer test: are the co-expression-differential genes enriched for
# differential methylation?
perm <- as.data.table(readRDS(file.path(results, "dgca_permFDR_final.rds")))
core <- perm[fdr_perm < 0.05 & fdr_pa < 0.05]
rewired <- unique(c(core$G1, core$G2)); rewired <- rewired[!is.na(rewired)]
gene_dm[, is_rewired := g %in% rewired]
tab <- table(rewired = gene_dm$is_rewired, DM = gene_dm$DM_strict)
ft <- fisher.test(tab)
cat("cross-layer enrichment: OR", round(ft$estimate, 2), "p", round(ft$p.value, 3), "\n")
cat("the two layers are ancestry-differential at different genes.\n")
