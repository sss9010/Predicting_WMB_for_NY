# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R/workflowr project for the **NY Winter Malting Barley (NY_WMB)** breeding program data analysis. It covers genomic and phenomic prediction of agronomic traits and malt quality in barley populations (DH = doubled haploids, RIL = recombinant inbred lines, Winter lines), using multi-year field trials and multi-spectral drone imagery.

## workflowr Commands

All analysis is managed via [workflowr](https://workflowr.github.io/workflowr/). Run these from the R console inside the project:

```r
# Build (render) a single analysis file
wflow_build("analysis/1. Data_Pre_Processing.Rmd")

# Build all analysis files
wflow_build()

# Publish (commit rendered HTML + source) to git
wflow_publish("analysis/1. Data_Pre_Processing.Rmd", "commit message")

# Check project status
wflow_status()
```

The `.Rprofile` automatically loads `workflowr` when the project is opened.

## Analysis Pipeline

The numbered `.Rmd` files in `analysis/` follow a sequential pipeline:

| File | Purpose |
|------|---------|
| `E1. Multi_spec_pre_processing.Rmd` | Import QGIS zonal statistics from drone flights, compute vegetation indices (NDVI, NDRE), QC and merge with plot IDs |
| `E2. Spectral_data_BLUP_new.Rmd` | Fit spatial mixed models on spectral/VI data to obtain BLUPs |
| `E3. Trait_spatial_mod_tester.Rmd` | Test spatial correction models for individual traits |
| `1. Data_Pre_Processing.Rmd` | Load pedigree, genotype (HapMap), and phenotype data; merge all datasets |
| `2. Trait heritability analysis.Rmd` | Heritability estimation for VI, malt quality, and agronomic traits via ASReml |
| `3. Spectral_data_BLUE.Rmd` | Compute BLUEs for spectral traits |
| `4. Exploratory_pheno_ analysis.Rmd` | Exploratory phenotypic analysis and visualization |
| `5. F_PCA.Rmd` | PCA on phenotypic/spectral data |
| `6. Phenomic kernel_build.Rmd` | Build phenomic relationship matrices from VI data for genomic prediction |
| `7. PP_ag_mq.Rmd` | Phenomic prediction for agronomic and malt quality traits |
| `8. GP_ag_mq.Rmd` | Genomic prediction (GBLUP) for agronomic and malt quality traits |
| `9. MT_GP.Rmd` | Multi-trait genomic prediction |
| `10. master_models.Rmd` | Combined/master prediction models |
| `11.DK_GP.Rmd` | Double-kernel genomic prediction |
| `12. Prediction_plots.Rmd` | Visualization of prediction results |

## Key R Packages

- **asreml** — mixed model fitting (spatial models, BLUPs/BLUEs, heritability)
- **ASRgenomics** — genomic relationship matrix construction, marker QC
- **ASRtriala** — trial analysis utilities
- **rrBLUP** — ridge regression BLUP for genomic prediction
- **workflowr** — reproducible research workflow and website publishing
- **tidyverse** (dplyr, tidyr, ggplot2, readr, readxl) — data wrangling and visualization

## Data Structure

```
data/
  geno/           # Genotype files
    *.hmp.txt     # HapMap format marker data (DH, RIL, parents, combined)
    *.vcf         # VCF format
    PLINK_*/      # PLINK .ped/.map files
    *.Rdata       # Pre-processed: Gmat, Relationship_matrix, wmb_GD_rrblup, etc.
  pheno/
    reg/          # Regular trial phenotypes (agronomic BLUPs/BLUEs, yield models)
    master/       # Master combined phenotype datasets
    MQ/           # Malt quality data
    VImstr/       # Vegetation index master datasets (NDVI, NDRE by year)
  Flights/        # QGIS zonal statistics CSVs from drone flights (y1–y4)
output/           # Processed/intermediate Rdata outputs
analysis/         # R Markdown analysis scripts (rendered to docs/)
docs/             # Rendered HTML (served as GitHub Pages site)
```

## Path Resolution Pattern

Scripts resolve data paths across environments using a fallback vector:

```r
base_path <- if (.Platform$OS.type == "windows") {
  "C:/Users/Siim Sepp/NY_WMB"
} else {
  "~/Documents/GitHub/NY_WMB"
}
```

Some newer scripts use a three-path fallback: Windows absolute, Unix absolute, and relative `"data"`.

## Genotype Data Notes

- Barley 50K SNP array data in HapMap format
- SNP positions mapped to **Morex V3** assembly (see `SNPPositions_50k_on_MorexV3.xlsx` and `barley50kMarkerPositions_Morex2019Assembly.gff`)
- Populations: DH (doubled haploids), RIL (recombinant inbred lines), Winter 2019 lines, and combined `all_filter.hmp.txt`
- Genomic relationship matrices pre-computed and stored as `Gmat.Rdata`, `Relationship_matrix.Rdata`, `wmb_GD_rrblup.Rdata`
