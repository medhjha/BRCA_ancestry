# Stage 1: download TCGA-BRCA data and set up the project.
# Pulls RNA, methylation, and copy number from the GDC, plus clinical data.
# The ancestry file (Carrot-Zhang 2020) is parsed separately in python.

library(TCGAbiolinks)
library(SummarizedExperiment)
library(yaml)

cfg <- read_yaml("config.yaml")
raw <- file.path(cfg$data, "raw", "tcga")
interim <- file.path(cfg$data, "interim")
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
dir.create(interim, recursive = TRUE, showWarnings = FALSE)

# RNA-seq (STAR counts). This is the slow one, about 1200 files.
q_rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(q_rna, directory = raw, method = "api", files.per.chunk = 20)
rna <- GDCprepare(q_rna, directory = raw)
saveRDS(rna, file.path(interim, "tcga_rna_se.rds"))

# Methylation (450K beta values). Large, about 900 samples.
q_meth <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)
GDCdownload(q_meth, directory = raw, method = "api", files.per.chunk = 10)
meth <- GDCprepare(q_meth, directory = raw)
saveRDS(meth, file.path(interim, "tcga_meth_se.rds"))

# Copy number. The product name matters here; "Gene Level Copy Number" is the
# current one, and there can be more than one workflow, so we pick ASCAT3.
q0 <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number"
)
workflows <- unique(getResults(q0)$analysis_workflow_type)
wf <- if ("ASCAT3" %in% workflows) "ASCAT3" else workflows[1]
q_cna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  workflow.type = wf
)
GDCdownload(q_cna, directory = raw, method = "api", files.per.chunk = 50)
cna <- GDCprepare(q_cna, directory = raw)
saveRDS(cna, file.path(interim, "tcga_cna_se.rds"))

cat("done. RNA, methylation, and copy number saved to", interim, "\n")
