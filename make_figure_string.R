# Figure 6: STRING enrichment and STRING confidence versus ancestry difference.
# Reads phyloframe_benchmark.csv written by 11_string_benchmark.R.

library(ggplot2)
library(data.table)
library(yaml)
library(patchwork)

cfg <- read_yaml("config.yaml")
results <- cfg$results
figdir <- file.path(results, "figures")
bm <- fread(file.path(results, "phyloframe_benchmark.csv"))

theme_paper <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 9.5, color = "grey40"),
        legend.position = "none")

# the two rates and the odds ratio come from the output of 11_string_benchmark.R
rates <- data.table(group = factor(c("Clean-core edges", "Random pairs\n(same genes)"),
                                   levels = c("Clean-core edges", "Random pairs\n(same genes)")),
                    pct = c(9.8, 1.3))
pa <- ggplot(rates, aes(group, pct, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(pct, "%")), vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("#2e6da4", "grey70")) +
  scale_y_continuous(limits = c(0, 12), expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Clean-core edges are present in STRING",
       subtitle = "odds ratio 8.4 compared with random pairs of the same genes",
       x = NULL, y = "percent present as a STRING edge") +
  theme_paper

insub <- bm[!is.na(string_score) & is.finite(zScoreDiff)]
insub[, absz := abs(zScoreDiff)]
pb <- ggplot(insub, aes(string_score, absz)) +
  geom_point(color = "#b03030", size = 1.3, alpha = 0.55) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30", linewidth = 0.6) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5, size = 3.8,
           label = "Spearman = 0.03\np = 0.60, n = 331") +
  labs(title = "STRING confidence does not predict the ancestry difference",
       subtitle = "clean-core edges present in STRING",
       x = "STRING combined confidence score", y = "absolute zScoreDiff") +
  theme_paper

fig <- pa + pb
ggsave(file.path(figdir, "phyloframe_orthogonality.png"), fig, width = 11, height = 5, dpi = 300, bg = "white")
ggsave(file.path(figdir, "phyloframe_orthogonality.pdf"), fig, width = 11, height = 5)
cat("saved Figure 6\n")
