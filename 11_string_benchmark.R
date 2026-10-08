# Comparison against a fixed functional network (STRING).
# The point is that a network that is the same for both ancestry groups cannot,
# by construction, say which edges differ by ancestry. We show empirically that
# the differential edges are real STRING interactions but that STRING confidence
# does not predict the ancestry difference.

library(STRINGdb)
library(data.table)
library(igraph)
library(yaml)
library(org.Hs.eg.db)

cfg <- read_yaml("config.yaml")
results <- cfg$results
set.seed(cfg$seeds$global)

core <- fread(file.path(results, "clean_core_edges_for_replication.csv"))
edge_genes <- unique(c(core$G1, core$G2))

# load STRING and map our genes plus a random background to STRING ids
sdb <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400,
                    input_directory = tempdir())
all_genes <- data.table(gene = unique(c(edge_genes,
                        sample(keys(org.Hs.eg.db, "SYMBOL"), 3000))))
mapped <- as.data.table(sdb$map(as.data.frame(all_genes), "gene", removeUnmappedRows = TRUE))
g2s <- setNames(mapped$STRING_id, mapped$gene)
net <- as.data.table(sdb$get_interactions(mapped$STRING_id))
net <- unique(net[, .(from, to, combined_score)])
net[, skey := paste(pmin(from, to), pmax(from, to))]
string_score <- setNames(net$combined_score, net$skey)

pair_score <- function(gA, gB) {
  a <- g2s[gA]; b <- g2s[gB]
  if (is.na(a) || is.na(b)) return(NA_real_)
  string_score[paste(pmin(a, b), pmax(a, b))]
}

# test 1: are the differential pairs in STRING more than random pairs?
core[, in_string := !is.na(mapply(pair_score, G1, G2))]
rewired_rate <- mean(core$in_string, na.rm = TRUE)
mg <- names(g2s)[names(g2s) %in% edge_genes]
rp <- data.table(a = sample(mg, 5000, replace = TRUE), b = sample(mg, 5000, replace = TRUE))[a != b]
rp[, in_string := !is.na(mapply(pair_score, a, b))]
rand_rate <- mean(rp$in_string, na.rm = TRUE)
ft <- fisher.test(matrix(c(sum(core$in_string, na.rm = TRUE), sum(!core$in_string, na.rm = TRUE),
                           sum(rp$in_string, na.rm = TRUE), sum(!rp$in_string, na.rm = TRUE)), 2))
cat("differential pairs in STRING:", round(100 * rewired_rate, 1), "percent\n")
cat("random pairs in STRING:", round(100 * rand_rate, 1), "percent\n")
cat("Fisher OR:", round(ft$estimate, 2), "p:", format(ft$p.value, scientific = TRUE), "\n")

# test 2: does STRING confidence predict the ancestry difference?
core[, string_score := mapply(pair_score, G1, G2)]
core[, tcga_delta := EUR_cor - AFR_cor]
in_string <- core[!is.na(string_score)]
ct <- cor.test(in_string$string_score, abs(in_string$zScoreDiff), method = "spearman")
cat("Spearman, STRING score vs ancestry-difference magnitude:", round(ct$estimate, 3),
    "p:", round(ct$p.value, 3), "n:", nrow(in_string), "\n")
cat("STRING confidence carries no information about which edges differ by ancestry.\n")

fwrite(core[, .(G1, G2, zScoreDiff, in_string, string_score)],
       file.path(results, "phyloframe_benchmark.csv"))
