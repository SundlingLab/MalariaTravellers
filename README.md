## R code repository for

# **Systems analysis shows a role of cytophilic antibodies in shaping innate tolerance to malaria**

[Lautenbach et al. Cell Reports 2022](https://doi.org/10.1016/j.celrep.2022.110709)

### A systems immunology analysis on natural malaria sheds light on disease tolerance mechanisms associated with cytophilic antibodies and gamma delta T cell expansion.

## Table of contents

-   [Repo description](#repo-description)
-   [Omics integration](#omics_integration)
-   [Figures](#figures)

## Repo description

This repository contains R code that was used for analysis and visualization.

This project used multiple omics data:

-   Plasma protein expression (Olink - NPX values)
-   Cell counts (Facs - White blood cell count normalized abundance)
-   IgG subclass response (Cumulative response score based on data from [Yman et al. 2019](https://doi.org/10.1186/s12916-019-1255-3))

## Omics Integration

### MEFISTO

MEFISTO is a Method for the Functional Integration of Spatial and Temporal Omics data and was used in this study in order to deconvolute the main sources of variation along the temporal axis of time after symptom onset for Plasma protein and cell abundance data mentioned above.

For more information, please read: [Velten et al. 2020](https://doi.org/10.1038/s41592-021-01343-9)

MEFISTO is part of the MultiOmicsFactorAnalysis (MOFA), publicly accessible here: [MEFISTO](https://biofam.github.io/MOFA2/MEFISTO.html)

-   `scripts/MEFISTO.R` Script used to train MEFISTO model

## Figures

R markdown scripts to reproduce Figures

-   `scripts/Figure1.Rmd` used for **Figure 1B, D-E**
-   `scripts/MEFISTO_downstream.Rmd` used for **Figures 2 & 4A,C & S2A-C & S4B**
-   `scripts/InternalExternalComparison.Rmd` used for **Figures S1C-D & S2D-G & S3**
-   `scripts/Correlation.Rmd` used for **Figures 3 & 4B & 7B-E**
-   `scripts/CRScore_calculation.Rmd` used for **Figure S7A-B**
-   `scripts/Additional_analysis.Rmd` used for **Figure S4A and S7C-H**
