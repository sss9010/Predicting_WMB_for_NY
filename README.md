# NY Winter Malting Barley — Genomic & Phenomic Prediction Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language: R](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)
[![Framework: workflowr](https://img.shields.io/badge/Framework-workflowr-brightgreen)](https://workflowr.github.io/workflowr/)

A reproducible R analysis pipeline for predicting **agronomic** and **malt quality** traits in the NY Winter Malting Barley (WMB) breeding program, combining UAV multispectral drone imagery, Barley 50K SNP genotype data, and multi-year field phenotypes across five environments (2021–2025).

More info: <https://cals.cornell.edu/field-crops/small-grains/small-grains-breeding-and-genetics>

---

> **TLDR:** Malt quality analysis is a considerable cost in barley breeding programs. Can UAVs be used to accelerate genetic gain and predict breeding values and phenotypes within field? Yes, but — phenomic prediction works better for agronomic traits and within years with more severe winter stress. Genomic prediction benefits from including correlated traits and varies in effectiveness depending on the trait. Many factors must be weighed before integrating these approaches into breeding programs, but the signal is there.

---

## Overview

This project develops and benchmarks a suite of prediction models for barley breeding selection, using three complementary data layers:

- **Spectral:** 5-band UAV multispectral imagery (blue, green, red, red-edge, NIR) at multiple timepoints per season
- **Genotype:** Barley 50K SNP array mapped to the Morex V3 assembly (DH, RIL, and Winter populations)
- **Phenotype:** Yield, heading date, and 12 malt quality traits (extract, protein, β-glucan, FAN, diastatic power, etc.)

Prediction approaches range from classical mixed models to machine learning and deep learning:

| Approach | Methods |
|----------|---------|
| Phenomic prediction | Linear regression, Ridge, LASSO, Elastic Net, PCR, Random Forest |
| Genomic prediction | GBLUP with spatial AR1×AR1 correction |
| Multi-trait genomic prediction | Bivariate GBLUP using spectral VIs as correlated helper traits |
| Double-kernel | Combined genomic + phenomic relationship kernels |
| DNN Synthetic VI | Bottleneck DNN learning optimal band combinations per trait |

---

## More Results with this data: DNN-Derived Synthetic Vegetation Index

A deep neural network with a **spectral bottleneck** (5 bands → 16 → 8 → **1 Synthetic VI** → trait) learns a data-optimised vegetation index for each target trait — outperforming fixed hand-crafted indices like NDVI under leave-one-environment-out cross-validation.

🔗 **Standalone DNN repository:** [sss9010/dnn-synthetic-vi](https://github.com/sss9010/dnn-synthetic-vi)  
📄 **Live report:** [DNN_Synthetic_VI.html](https://sss9010.github.io/dnn-synthetic-vi/DNN_Synthetic_VI.html)

---

## Analysis Pipeline

| Script | Description |
|--------|-------------|
| `1. Trait heritability analysis.Rmd` | Estimate h² and H² for VIs, malt quality, and agronomic traits via ASReml |
| `2. Exploratory_pheno_ analysis.Rmd` | Exploratory phenotypic analysis and visualisation |
| `3. Spectral_data_BLUE.Rmd` | AR1×AR1 spatial BLUE models for spectral and agronomic traits; build VI time-series matrix |
| `4. F_PCA.Rmd` | PCA and FPCA on spectral time-series; extract functional principal components |
| `5. Phenomic kernel_build.Rmd` | Build phenomic relationship matrices from VI time-series data |
| `6. PP_ag_mq.Rmd` | Phenomic prediction for agronomic and malt quality traits |
| `7. GP_ag_mq.Rmd` | GBLUP cross-validation for agronomic and malt quality traits |
| `8. MT_GP.Rmd` | Multi-trait GBLUP using spectral VIs and FPCs as correlated helper traits |
| `10.DK_GP.Rmd` | Double-kernel GBLUP (genomic G + phenomic P kernels) |
| `11. Prediction_plots.Rmd` | Visualisation and cross-method comparison of prediction accuracies |

---

## Data

| Location | Contents |
|----------|----------|
| `data/geno/` | 50K SNP array in HapMap, VCF, and PLINK formats (Morex V3); DH, RIL, Winter populations |
| `data/Flights/` | QGIS zonal statistics CSVs from UAV flights across five environments |
| `data/ag_BLUE_spatial.RData` | Agronomic BLUEs (spatially corrected via AR1×AR1) |
| `data/spec_BLUE_spatial.RData` | Spectral BLUEs across timepoints |
| `data/data_mq.Rdata` | Malt quality phenotype data |
| `data/MT_pred_mat.Rdata` / `MT_raw_mat.Rdata` | Phenomic relationship matrices (BLUE & raw VI) |
| `data/F_PCA_scoresG.Rdata` | Functional PCA scores per VI and environment |
| `data/WMB_pheno.Rdata` | Master merged dataset (bands + VIs + traits) |

**Environments:** HELF24, KET21, MCG23, MCG25, SNY22  
**Populations:** DH (doubled haploids), RIL (recombinant inbred lines)

> `.Rdata`  input files are available in the data folder.

---

## Dependencies

```r
# Mixed models & genomics (asreml requires a VSNi licence — free for academia)
install.packages(c("asreml", "ASRgenomics", "ASRtriala", "rrBLUP"))

# Machine learning
install.packages(c("glmnet", "ranger", "caret"))

# Functional data analysis
install.packages("fdapace")

# Deep learning (DNN Synthetic VI)
install.packages("torch"); torch::install_torch()

# Data wrangling & visualisation
install.packages(c("tidyverse", "ggplot2", "patchwork", "corrplot",
                   "readxl", "writexl", "knitr", "kableExtra"))

# Reproducibility
install.packages("workflowr")
```

---

## Reproducing the Analysis

All analysis is managed via [workflowr](https://workflowr.github.io/workflowr/). From an R console in the project root:

```r
wflow_build("analysis/7. GP_ag_mq.Rmd")  # render a single step
wflow_build()                              # render full pipeline
wflow_publish("analysis/7. GP_ag_mq.Rmd", "message")  # publish to docs/
wflow_status()                             # check project status
```

---

## Citation

```
Sepp, S.S., Sorrells, M.E. (2026). NY Winter Malting Barley —
Genomic & Phenomic Prediction Pipeline. GitHub.
https://github.com/sss9010/Predicting_WMB_for_NY
```

A machine-readable citation is in [`CITATION.cff`](CITATION.cff).

---

## Contact

**Siim Sepp** — sss322@cornell.edu  
NY Winter Malting Barley Breeding Programme, Cornell University
