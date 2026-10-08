# Export the clean-core edge list that several downstream scripts read.
# Run this once after 05b (permutation FDR) and 07 (composition control).

library(data.table)
library(yaml)

cfg <- read_yaml("config.yaml")
results <- cfg$results

perm <- as.data.table(readRDS(file.path(results, "dgca_permFDR_final.rds")))
core <- perm[fdr_perm < 0.05 & fdr_pa < 0.05]

# join the correlations from the parametric result
full <- as.data.table(readRDS(file.path(results, "dgca_final_EUR_vs_AFR.rds")))
full[, key := paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2))]
out <- merge(core[, .(key, fdr_perm)],
             full[, .(key, G1, G2, EUR_cor, AFR_cor, zScoreDiff)],
             by = "key")
out <- out[!is.na(G1) & !is.na(G2)]

fwrite(out, file.path(results, "clean_core_edges_for_replication.csv"))
cat("exported", nrow(out), "clean-core edges\n")
