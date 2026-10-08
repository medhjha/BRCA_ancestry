# Figure 3: hub network of clean-core edges.
# Reads the clean-core edge export (06_export_clean_core.R). Draws the top hub
# genes, with node size by edge count and edge color by direction (red tighter in
# AFR, blue tighter in EUR).

library(data.table)
library(igraph)
library(yaml)

cfg <- read_yaml("config.yaml")
results <- cfg$results
figdir <- file.path(results, "figures")
dir.create(figdir, showWarnings = FALSE)

core <- fread(file.path(results, "clean_core_edges_for_replication.csv"))
core[, afr_tighter := abs(AFR_cor) > abs(EUR_cor)]

# keep edges among the top hub genes so the picture is readable
deg <- rbind(core[, .(g = G1)], core[, .(g = G2)])[, .N, by = g][order(-N)]
top_hubs <- deg[1:40, g]
sub <- core[G1 %in% top_hubs & G2 %in% top_hubs]

g <- graph_from_data_frame(
  sub[, .(from = G1, to = G2, afr_tighter)], directed = FALSE)
deg_map <- setNames(deg$N, deg$g)
V(g)$size <- 4 + 2 * sqrt(deg_map[V(g)$name])
E(g)$color <- ifelse(sub$afr_tighter, "#c0392b", "#2e6da4")

png(file.path(figdir, "rewired_network.png"), width = 2400, height = 2000, res = 300)
set.seed(1)
plot(g, layout = layout_with_fr(g),
     vertex.label.cex = 0.6, vertex.label.color = "black",
     vertex.color = "grey90", vertex.frame.color = "grey50",
     edge.width = 1.2, main = "Ancestry-differential co-expression hubs")
legend("bottomleft", legend = c("tighter in AFR", "tighter in EUR"),
       col = c("#c0392b", "#2e6da4"), lwd = 2, bty = "n", cex = 0.7)
dev.off()
cat("saved Figure 3 (rewired_network.png)\n")
