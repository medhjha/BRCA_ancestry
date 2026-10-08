# Run order

Run the scripts in this order. Each one reads what the previous ones saved, so
they can be run one at a time. Set the paths in config.yaml first.

## Setup

Edit config.yaml so the root, data, and results paths point at your machine.
The R scripts expect config.yaml in the working directory.

## Data

    python/parse_ancestry.py        download and parse the ancestry table
    R/01_download_tcga.R            download RNA, methylation, copy number
    R/02_cohort_audit.R             freeze the 3-layer cohort

## Preprocessing

    R/03_preprocess.R               VST, M-values, tumor-half copy number
    R/03b_covariates_and_batch.R    covariates, and the batch confounding check

## Network and differential co-expression

    R/04_wgcna_and_preservation.R   WGCNA on EUR, module preservation to AFR
    R/05_differential_coexpression.R  clean gene universe, DGCA, parametric FDR
    R/07_composition_control.R      composition control, clean-core, enrichment
    R/05b_permutation_fdr.R         permutation FDR (needs 07 for the clean-core)
    R/06_export_clean_core.R        export the clean-core edge list
    R/12_directional_enrichment.R   hub genes by direction, GO enrichment

Note: 07 writes the composition classification that 05b reads, so run 07 before
05b even though the numbering suggests otherwise.

## Confound controls and robustness

    R/08_mixedsite_robustness.R     re-run within mixed sites
    R/09_confound_controls.R        subtype, age, and stage controls

## Second layer and validation

    R/10_methylation.R              differential methylation, cross-layer test
    python/cptac_replication.py     CPTAC RNA replication and protein
    R/11_string_benchmark.R         comparison against STRING

## Figures

    figures/make_figures.R          main figures
    figures/make_supplementary_figures.R  supplementary figures S1 to S4
    figures/make_figure_string.R    Figure 6 (STRING)
    figures/make_figure_network.R   Figure 3 (hub network)


## A note on the clean-core edge export

Several scripts read a file called clean_core_edges_for_replication.csv, which is
the clean-core edge list (columns G1, G2, EUR_cor, AFR_cor, zScoreDiff, fdr_perm).
Export it once after 05b and 07 with a short script that filters
dgca_permFDR_final.rds to fdr_perm < 0.05 and fdr_pa < 0.05 and joins the
correlations from dgca_final_EUR_vs_AFR.rds.
