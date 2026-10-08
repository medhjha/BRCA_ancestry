# Make the main and supplementary figures from the saved result files.
# Each figure reads one or two result files and writes a png and a pdf.

library(ggplot2)
library(data.table)
library(yaml)
library(patchwork)

cfg <- read_yaml("config.yaml")
results <- cfg$results
figdir <- file.path(results, "figures")
dir.create(figdir, showWarnings = FALSE)

theme_paper <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 9.5, color = "grey40"),
        plot.caption = element_text(size = 8, color = "grey50", hjust = 0))

save_fig <- function(plot, name, w = 7, h = 5.5) {
  ggsave(file.path(figdir, paste0(name, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(figdir, paste0(name, ".pdf")), plot, width = w, height = h)
}

# Figure 1: module preservation
pres <- fread(file.path(results, "module_preservation_table.csv"))
pres[, module := factor(module)]
f1 <- ggplot(pres, aes(moduleSize, Zsummary)) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "orange") +
  geom_point(aes(color = module), size = 4) +
  labs(title = "EUR co-expression modules are preserved in AFR",
       subtitle = "Module-level architecture is ancestry-invariant",
       x = "module size (genes)", y = "preservation Zsummary (EUR to AFR)",
       caption = "All modules above Zsummary 10 are strongly preserved.") +
  theme_paper + theme(legend.position = "none")
save_fig(f1, "module_preservation")

# Figure 2: DGCA volcano
dc <- as.data.table(readRDS(file.path(results, "dgca_final_EUR_vs_AFR.rds")))
dc[, neglogFDR := -log10(pmax(fdr, 1e-300))]
dc[, sig := fdr < 0.05]
f2 <- ggplot(dc, aes(zScoreDiff, neglogFDR)) +
  geom_point(aes(color = Classes, alpha = sig), size = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_alpha_manual(values = c(`TRUE` = 0.85, `FALSE` = 0.12), guide = "none") +
  labs(title = "Ancestry-differential co-expression edges",
       x = "zScoreDiff (positive: higher correlation in AFR, negative: higher in EUR)",
       y = "-log10(FDR)", color = "class") +
  theme_paper
save_fig(f2, "dgca_volcano", w = 8, h = 6)

# Figure 3: mixed-site concordance
cmp <- fread(file.path(results, "dgca_full_vs_mixedsite.csv"))
cmp[, sig := fdr_full < 0.05]
sp <- cor(cmp$z_full, cmp$z_mixed, method = "spearman", use = "complete.obs")
f3 <- ggplot(cmp, aes(z_full, z_mixed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = sig), size = 0.5, alpha = 0.4) +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70"),
                     labels = c("ns", "FDR<0.05"), name = NULL) +
  annotate("text", x = min(cmp$z_full, na.rm = TRUE) * 0.7,
           y = max(cmp$z_mixed, na.rm = TRUE) * 0.85,
           label = sprintf("Spearman = %.3f", sp), hjust = 0, size = 4, fontface = "bold") +
  labs(title = "Differential signal survives site control",
       x = "zScoreDiff (full)", y = "zScoreDiff (site-controlled)") +
  theme_paper
save_fig(f3, "mixedsite_concordance", w = 6.5, h = 6)

# Figure 5: CPTAC replication
rep_rna <- fread(file.path(results, "cptac_replication.csv"))
rep_rna[, td := EUR_cor - AFR_cor]
rep_rna[, concordant := sign(td) == sign(cpd)]
f5 <- ggplot(rep_rna, aes(td, cpd)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(aes(color = concordant), size = 0.6, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "firebrick") +
  scale_color_manual(values = c(`TRUE` = "#1e7d34", `FALSE` = "grey65")) +
  labs(title = "RNA replication: TCGA to CPTAC",
       x = "TCGA delta-corr (EUR minus AFR)", y = "CPTAC delta-corr (EUR minus AFR)") +
  theme_paper + theme(legend.position = "none")
save_fig(f5, "replication_scatter", w = 6.5, h = 5.5)

cat("figures saved to", figdir, "\n")
