# Ancestry-differential co-expression in breast cancer

Code for the analysis of ancestry-differential gene co-expression networks in
breast cancer, comparing European-ancestry (EUR) and African-ancestry (AFR)
tumors in TCGA, with replication in CPTAC.

## What the analysis does

We test whether gene co-expression in breast cancer differs between EUR and AFR
patients, and at what scale. The short version of the result: the large
co-expression modules are the same in both groups, but individual gene-gene
edges differ, and those differences hold up under composition, subtype, age,
stage, and site controls, replicate in an independent cohort, and are not
captured by a fixed functional network.

## Layout

    R/          R scripts for each analysis stage, numbered in run order
    python/     Python scripts (CPTAC data, ancestry parsing, replication)
    figures/    R scripts that make the main and supplementary figures
    config.yaml example config; every script reads paths and parameters from here

## How it is organized

The pipeline runs in stages. Each script reads its inputs from disk (saved by an
earlier stage) and writes its outputs back to disk, so stages can be run one at a
time. Paths and parameters live in config.yaml; change the paths there to match
where your data lives.

    Stage 1   download TCGA data and parse ancestry calls
    Stage 2   audit which patients have all layers, freeze the cohort
    Stage 3   preprocess expression, methylation, copy number, and covariates
    Stage 4   build the co-expression network (WGCNA) on EUR tumors
    Stage 5   module preservation, differential co-expression (DGCA), gradient
    Stage 7   biological interpretation and composition control
    Stage 8   mixed-site robustness
    replication   CPTAC RNA replication and protein corroboration
    methylation   differential methylation and the cross-layer test
    benchmark     comparison against a fixed functional network (STRING)

## Requirements

R 4.4 with Bioconductor 3.20. Main packages: TCGAbiolinks, WGCNA, DGCA, DESeq2,
limma, sesame, genefu, STRINGdb, lionessR, clusterProfiler, estimate,
immunedeconv. Python 3.11 with cptac, pandas, numpy, scipy.

Exact versions used are listed in the manuscript methods section.

## Repository URL

Replace the placeholder repository URL in the manuscript methods (software and
data availability) with this repository address once it is published.

## Data

All input data is public. TCGA breast cancer data comes from the NCI Genomic
Data Commons through TCGAbiolinks. CPTAC breast cancer data comes from the cptac
Python package. Ancestry calls come from Carrot-Zhang et al. 2020.

## Note on the filesystem

The scripts were written to run off an external drive, so paths point at a mount
point set in config.yaml. On some filesystems (for example exFAT) conda
environments and R package libraries do not work reliably, so keep those on a
local disk and keep only data and results on the external drive if you use one.
