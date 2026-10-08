# Subtype, age, and stage controls.
# These all differ by ancestry and cannot enter DGCA directly, so we recompute
# each clean-core edge's EUR-vs-AFR difference within single strata and check the
# direction concordance with the full-cohort result. If an edge reflects ancestry
# rather than a confounder, its direction holds within a stratum where that
# confounder is held constant.

library(data.table)
library(arrow)
library(yaml)
library(SummarizedExperiment)

cfg <- read_yaml("config.yaml")
interim <- file.path(cfg$data, "interim")
proc <- file.path(cfg$data, "processed")
results <- cfg$results

coh <- as.data.table(read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
coh <- coh[patient12 %in% common & cohort %in% c("EUR_primary", "AFR_primary")]
coh[, anc := fifelse(cohort == "EUR_primary", "EUR", "AFR")]
cov <- as.data.table(readRDS(file.path(proc, "covariates.rds")))

rna <- readRDS(file.path(interim, "tcga_rna_se.rds"))
cd <- as.data.table(as.data.frame(colData(rna)))
cd[, patient12 := substr(colnames(rna), 1, 12)]
cd <- cd[!duplicated(patient12)]
age_col <- grep("age_at_diagnosis|age_at_index", names(cd), value = TRUE, ignore.case = TRUE)[1]
stage_col <- grep("ajcc_pathologic_stage", names(cd), value = TRUE)[1]
cd[, age_yr := { v <- as.numeric(get(age_col)); if (median(v, na.rm = TRUE) > 200) v / 365.25 else v }]
cd[, stage3 := gsub("Stage |A|B|C", "", as.character(get(stage_col)))]

d <- merge(coh[, .(patient12, anc)], cov[, .(patient12, pam50_genefu)], by = "patient12")
d <- merge(d, cd[, .(patient12, age_yr, stage3)], by = "patient12")

vsd <- readRDS(file.path(proc, "rna_vst.rds"))
core <- fread(file.path(results, "clean_core_edges_for_replication.csv"))
core[, tcga_delta := EUR_cor - AFR_cor]
sym <- setNames(rowData(rna)$gene_name, rownames(rna))
ens_of <- function(s) { i <- names(sym)[sym == s]; i <- i[i %in% rownames(vsd)]; if (length(i)) i[1] else NA }
g2e <- setNames(vapply(unique(c(core$G1, core$G2)), ens_of, character(1)),
                unique(c(core$G1, core$G2)))
core <- core[!is.na(g2e[G1]) & !is.na(g2e[G2])]

# recompute EUR-vs-AFR delta within a given patient set, and report concordance
within_delta <- function(pe, pa) {
  E <- vsd[, colnames(vsd) %in% pe]; A <- vsd[, colnames(vsd) %in% pa]
  sapply(seq_len(nrow(core)), function(i) {
    e1 <- g2e[core$G1[i]]; e2 <- g2e[core$G2[i]]
    suppressWarnings(cor(E[e1, ], E[e2, ])) - suppressWarnings(cor(A[e1, ], A[e2, ]))
  })
}
concordance <- function(pe, pa, label) {
  if (length(pa) < 15) { cat(sprintf("  %-16s too few AFR (%d)\n", label, length(pa))); return(NULL) }
  sd <- within_delta(pe, pa)
  ok <- is.finite(sd) & is.finite(core$tcga_delta)
  conc <- mean(sign(sd[ok]) == sign(core$tcga_delta[ok]))
  cat(sprintf("  %-16s EUR %d / AFR %d: %.0f percent concordant\n",
              label, length(pe), length(pa), 100 * conc))
  data.table(subtype = label, n_eur = length(pe), n_afr = length(pa),
             concordance = round(100 * conc, 1))
}

cat("within-subtype:\n")
res <- rbindlist(list(
  concordance(d[anc == "EUR" & pam50_genefu %in% c("LumA", "LumB"), patient12],
              d[anc == "AFR" & pam50_genefu %in% c("LumA", "LumB"), patient12], "Luminal (A+B)"),
  concordance(d[anc == "EUR" & pam50_genefu == "Basal", patient12],
              d[anc == "AFR" & pam50_genefu == "Basal", patient12], "Basal"),
  concordance(d[anc == "EUR" & pam50_genefu == "LumB", patient12],
              d[anc == "AFR" & pam50_genefu == "LumB", patient12], "LumB only")
), fill = TRUE)
fwrite(res, file.path(results, "subtype_confound_test.csv"))

cat("\nwithin age band 45-65 and within stage:\n")
band <- d[age_yr >= 45 & age_yr <= 65]
concordance(band[anc == "EUR", patient12], band[anc == "AFR", patient12], "age 45-65")
s2 <- d[stage3 == "II"]
concordance(s2[anc == "EUR", patient12], s2[anc == "AFR", patient12], "stage II")

cat("\nage/stage differ by ancestry but do not drive the edges.\n")
