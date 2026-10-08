# Stage 4 and 5a: build the co-expression network on EUR tumors and test whether
# the modules are preserved in AFR tumors.

library(WGCNA)
library(data.table)
library(arrow)
library(yaml)
library(matrixStats)

enableWGCNAThreads(4)
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

# pick the 5000 most variable genes in EUR, so the reference network is defined
# on the discovery group only
vsd_eur <- vsd[, colnames(vsd) %in% eur]
top <- order(rowVars(vsd_eur), decreasing = TRUE)[1:5000]
gene_set <- rownames(vsd_eur)[top]
dat_eur <- t(vsd_eur[top, ])

# choose the soft threshold power by scale-free fit
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(dat_eur, powerVector = powers,
                         networkType = "signed", corFnc = "bicor")
power <- sft$powerEstimate
if (is.na(power)) power <- 12
cat("chosen soft power:", power, "\n")

net <- blockwiseModules(dat_eur, power = power, networkType = "signed",
                        TOMType = "signed", corType = "bicor",
                        maxPOutliers = 0.10, minModuleSize = 30,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        maxBlockSize = 5000, verbose = 1)
cat("modules found:\n")
print(table(net$colors))

saveRDS(net, file.path(proc, "wgcna_net_eur.rds"))
saveRDS(gene_set, file.path(proc, "wgcna_genes.rds"))
saveRDS(split(gene_set, net$colors), file.path(proc, "module_membership.rds"))

# module preservation: do the EUR modules hold in AFR?
dat_afr <- t(vsd[gene_set, colnames(vsd) %in% afr])
multi_expr <- list(EUR = list(data = dat_eur), AFR = list(data = dat_afr))
multi_color <- list(EUR = net$colors)
mp <- modulePreservation(multi_expr, multi_color,
                         referenceNetworks = 1, testNetworks = 2,
                         nPermutations = cfg$params$perm_preservation,
                         randomSeed = cfg$seeds$global, verbose = 2)
saveRDS(mp, file.path(results, "module_preservation_EURtoAFR.rds"))

# the Zsummary table lives under a version specific path, so find the AFR test
# frame by name rather than by position
ref <- names(mp$preservation$Z)[1]
z_children <- names(mp$preservation$Z[[ref]])
test <- z_children[grepl("AFR", z_children) &
                     sapply(mp$preservation$Z[[ref]], is.data.frame)]
zdf <- mp$preservation$Z[[ref]][[test]]
obs <- mp$preservation$observed[[ref]][[test]]
res <- data.table(module = rownames(zdf), moduleSize = zdf$moduleSize,
                  Zsummary = zdf$Zsummary.pres, medianRank = obs$medianRank.pres)
res <- res[!(module %in% c("0", "0.1", "gold"))]
cat("\nmodule preservation (Zsummary > 10 is strong):\n")
print(res[order(Zsummary)])
fwrite(res, file.path(results, "module_preservation_table.csv"))
