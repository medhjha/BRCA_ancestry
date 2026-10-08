# Make the supplementary figures from saved result files.
# Numbering follows the order in which the figures are first cited in the text:
#   S1 batch confounding, S2 subtype control, S3 methylation, S4 cross-layer.

library(ggplot2)
library(data.table)
library(yaml)
library(patchwork)

cfg <- read_yaml("config.yaml")
results <- cfg$results
proc <- file.path(cfg$data, "processed")
figdir <- file.path(results, "figures")
dir.create(figdir, showWarnings = FALSE)

theme_paper <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 9.5, color = "grey40"),
        plot.caption = element_text(size = 8, color = "grey50", hjust = 0))

save_fig <- function(plot, name, w, h) {
  ggsave(file.path(figdir, paste0(name, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(figdir, paste0(name, ".pdf")), plot, width = w, height = h)
}

# S1: tissue source site and ancestry
meta <- as.data.table(readRDS(file.path(proc, "batch_ancestry_meta.rds")))
common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
meta <- meta[patient12 %in% common & anc2 %in% c("EUR", "AFR")]
sc <- meta[, .(EUR = sum(anc2 == "EUR"), AFR = sum(anc2 == "AFR")), by = tss][order(-(EUR + AFR))]
sc[, tss := factor(tss, levels = tss)]
scm <- melt(sc, id.vars = "tss", variable.name = "ancestry", value.name = "n")
s1 <- ggplot(scm, aes(tss, n, fill = ancestry)) +
  geom_col(position = "stack", width = 0.85) +
  scale_fill_manual(values = c(EUR = "#2e6da4", AFR = "#c0392b")) +
  labs(title = "Tissue collection site is confounded with ancestry",
       subtitle = "Cramer's V = 0.579 (site) and 0.585 (plate)",
       x = "tissue source site (ordered by size)", y = "patients", fill = NULL,
       caption = "Some sites are single-ancestry, so a site-restricted analysis is used instead of batch correction.") +
  theme_paper + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
save_fig(s1, "supp_S1_batch_confounding", 9, 4.5)

# S2: subtype control
sub <- fread(file.path(results, "subtype_confound_test.csv"))
sub[, lab := paste0(subtype, "\n(EUR ", n_eur, " / AFR ", n_afr, ")")]
sub[, lab := factor(lab, levels = lab[order(-concordance)])]
s2 <- ggplot(sub, aes(lab, concordance)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50") +
  geom_col(fill = "#238b45", width = 0.6) +
  geom_text(aes(label = paste0(concordance, "%")), vjust = -0.4, size = 4.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Differential edges keep their direction within molecular subtype",
       subtitle = "percent of edges whose EUR-vs-AFR direction matches the full cohort",
       x = NULL, y = "percent edges concordant",
       caption = "Dashed line = 50 percent (chance). Her2-enriched (underpowered) and Normal-like are not shown.") +
  theme_paper
save_fig(s2, "supp_S2_subtype_control", 7, 5)

# S3: methylation effect sizes and per-chromosome direction
tt <- as.data.table(readRDS(file.path(results, "dmp_probe_EUR_vs_AFR.rds")))
sig <- tt[fdr < 0.05]
s3a <- ggplot(sig, aes(abs(delta_M_AFRvsEUR))) +
  geom_histogram(bins = 50, fill = "#6a51a3", color = "white", linewidth = 0.1) +
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "grey30") +
  labs(title = "Effect sizes of differentially methylated probes",
       subtitle = "Significant promoter probes (FDR < 0.05); dashed line marks 0.2",
       x = "absolute delta M (AFR vs EUR)", y = "probes") +
  theme_paper
p2g <- as.data.table(readRDS(file.path(proc, "meth_promoter_probe2gene.rds")))
sig2 <- merge(unique(sig[, .(probe, delta_M_AFRvsEUR)]), unique(p2g[, .(probe, chr)]), by = "probe")
chrsk <- sig2[, .(pct_hyper = 100 * mean(delta_M_AFRvsEUR > 0)), by = chr]
chrsk <- chrsk[grepl("^chr[0-9XY]+$", chr)]
chrsk[, ord := suppressWarnings(as.integer(gsub("chr|X|Y", "", chr)))]
chrsk <- chrsk[order(ord)]
chrsk[, chr := factor(chr, levels = chr)]
s3b <- ggplot(chrsk, aes(chr, pct_hyper)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50") +
  geom_col(fill = "#41ab5d", width = 0.8) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Direction of methylation differences by chromosome",
       subtitle = "percent of significant probes with higher methylation in AFR",
       x = NULL, y = "percent higher in AFR") +
  theme_paper + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))
save_fig(s3a / s3b, "supp_S3_methylation_effectsize", 8, 8)

# S4: cross-layer comparison
gs <- as.data.table(readRDS(file.path(results, "dmp_gene_strict.rds")))
perm <- as.data.table(readRDS(file.path(results, "dgca_permFDR_final.rds")))
core <- perm[fdr_perm < 0.05 & fdr_pa < 0.05]
rewired <- unique(c(core$G1, core$G2))
rewired <- rewired[!is.na(rewired)]
gs[, is_rewired := g %in% rewired]
dm <- gs[, .(dm_rate = 100 * mean(DM_strict)), by = is_rewired]
dm[, group := fifelse(is_rewired, "Co-expression\ndifferential genes", "Other\ngenes")]
s4 <- ggplot(dm, aes(group, dm_rate, fill = group)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = paste0(round(dm_rate), "%")), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Co-expression\ndifferential genes" = "#2e6da4", "Other\ngenes" = "grey70")) +
  scale_y_continuous(limits = c(0, 70), expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Co-expression and methylation differences are at different genes",
       subtitle = "percent differentially methylated (FDR < 0.05 and |delta M| >= 0.2); odds ratio 1.00",
       x = NULL, y = "percent differentially methylated") +
  theme_paper + theme(legend.position = "none")
save_fig(s4, "supp_S4_crosslayer_null", 7, 5)

cat("supplementary figures saved to", figdir, "\n")
