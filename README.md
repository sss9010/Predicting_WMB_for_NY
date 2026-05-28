# NY Winter Malting Barley — Genomic & Phenomic Prediction Pipeline

[![workflowr](https://img.shields.io/badge/workflowr-reproducible-brightgreen)](https://workflowr.github.io/workflowr/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A reproducible R analysis pipeline for predicting **agronomic** and **malt quality** traits in the NY Winter Malting Barley (WMB) breeding program, combining UAV multispectral drone imagery, Barley 50K SNP genotype data, and multi-year field phenotypes across five environments (2021–2025).

TLDR: Malt quality analysis is a considerable cost in barley breeding programs. Can UAVs be used to accelerate genetic gain and predict breeding values and phenotypes within field? Yes, but... Phenomic prediction works better for agronomic traits and within years with more severe winter stress. Genomic prediction benefits from including a correlated traits and varies in effectiveness depending on the trait. Many factors must be weighed before integrating these approaches into breeding programs, but the signal is there. Sometimes more, sometimes less. 

More info: https://cals.cornell.edu/field-crops/small-grains/small-grains-breeding-and-genetics
---

## Overview

This project develops and benchmarks a suite of prediction models for barley breeding selection, using three complementary data layers:

- **Spectral:** NDVI and PSRI from UAV multispectral flights at multiple timepoints per season
- **Genotype:** Barley 50K SNP array mapped to the Morex V3 assembly
- **Phenotype:** Yield, heading date, and 12 malt quality traits (extract, protein, β-glucan, FAN, diastatic power, etc.)

Prediction approaches range from classical mixed models to machine learning and deep learning:

| Approach | Methods |
|----------|---------|
| Phenomic prediction | Linear regression, Ridge, LASSO, Elastic Net, PCR, Random Forest |
| Genomic prediction | GBLUP with spatial AR1×AR1 correction |
| Multi-trait genomic prediction | Bivariate GBLUP using spectral VIs as correlated helper traits |
| Double-kernel | Combined genomic + phenomic relationship kernels |

---

## Analysis Pipeline

| File | Description |
|------|-------------|
| `1. Trait heritability analysis.Rmd` | Estimate h² and H² for VIs, malt quality, and agronomic traits via ASReml |
| `2. Exploratory_pheno_ analysis.Rmd` | Trait distributions, correlation matrices, heading date and flight timeline visualisation |
| `3. Spectral_data_BLUE.Rmd` | Spatial BLUE models for spectral and agronomic traits; build VI time-series matrix |
| `4. F_PCA.Rmd` | PCA and FPCA on spectral time-series; extract functional principal components |
| `5. Phenomic kernel_build.Rmd` | Build phenomic relationship matrices from VI time-series data |
| `6. PP_ag_mq.Rmd` | Phenomic prediction for agronomic (BLUE) and malt quality (raw) traits |
| `7. GP_ag_mq.Rmd` | GBLUP cross-validation for agronomic (BLUEs) and malt quality (raw) traits |
| `8. MT_GP.Rmd` | Multi-trait GBLUP using spectral VIs and FPCs as correlated helper traits |
| `10.DK_GP.Rmd` | Double-kernel GBLUP (genomic G + phenomic P kernels) |
| `11. Prediction_plots.Rmd` | Visualisation and cross-method comparison of prediction accuracies |

---

## Rendered Analyses

Full rendered HTML outputs are available at the project GitHub Pages site:
**https://sss9010.github.io/Predicting_WMB_for_NY**

---

## Dependencies

Core packages (install from CRAN unless noted):

```r
# Mixed models & genomics (asreml requires a VSNi licence — free for academia)
install.packages(c("asreml", "ASRgenomics", "ASRtriala", "rrBLUP"))

# Machine learning
install.packages(c("caret", "glmnet", "randomForest"))

# Functional data analysis
install.packages("fda")

# Data wrangling & visualisation
install.packages(c("tidyverse", "readxl", "writexl", "data.table", "purrr"))

# Reproducibility
install.packages("workflowr")
```

---

## Data

| Directory | Contents |
|-----------|----------|
| `data/geno/` | 50K SNP array in HapMap, VCF, and PLINK formats (Morex V3); DH, RIL, Winter populations |
| `data/pheno/` | Agronomic BLUPs/BLUEs, malt quality data, VI master datasets by year |
| `data/Flights/` | QGIS zonal statistics CSVs from UAV flights across five environments |

> Raw genotype files are not redistributed in this repository. Processed intermediate `.Rdata` files are available on request.

---

## Reproducing the Analysis

All analysis is managed via [workflowr](https://workflowr.github.io/workflowr/). From an R console in the project root:

```r
# Render a single step
wflow_build("analysis/1. Data_Pre_Processing.Rmd")

# Render the full pipeline
wflow_build()

# Publish rendered HTML to docs/ and commit
wflow_publish("analysis/1. Data_Pre_Processing.Rmd", "message")

# Check project status
wflow_status()
```

---

## Citation

If you use this pipeline or results, please cite:

> Sepp, S. (2025). *NY Winter Malting Barley Genomic and Phenomic Prediction Pipeline*. GitHub. https://github.com/sss9010/Predicting_WMB_for_NY

---

## Contact

Siim Sepp · sss@322cornell.edu
