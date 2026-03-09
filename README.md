# NY_WMB

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr


# 🌾 Winter Barley UAV Phenomics and Trait Prediction Pipeline

This repository documents the analysis pipeline for spectral, agronomic (AG), and malting quality (MQ) trait modeling in winter barley, based on UAV multispectral data, genotype, and field phenotypes.

---

## Overview of Analysis Steps

### **1. Data Pre-Processing for Spectral, Agronomic, and MQ Traits**

- Prepare spectral (VI), agronomic (AG), and malting quality (MQ) datasets into unified dataframes.
- Process genotype data and integrate pedigree information.

**Outputs:**
- Dataframes joinable by `plot_id` or `GID`.
  - `data/WMB_pheno_all.Rdata` — individual datasets
  - `data/WMB_pheno.Rdata` — merged dataset
- Vegetation index matrix: `data/VImat.Rdata`

---

s






### **2. Trait Heritability Analysis**

- Estimate:
  - Narrow-sense heritability (h²)
  - Broad-sense heritability (H²)
- Applies to:
  - Agronomic traits
  - Malting quality traits
  - Spectral indices (VIs)
  - PCA / FPCA traits

**Outputs:**
- `analysis/heritabilities.xlsx`
- Heritability bar plots

---








### **3. Spectral First-Stage BLUP Modeling Over Timepoints**

- Fit spatial mixed models with autoregressive terms for each VI at each timepoint.
- Extract BLUEs (Best Linear Unbiased Estimates) and BLUPs for downstream analysis.

**Outputs:**
- **Spectral BLUEs:**
  - `data/spec_BLUE_spatial_fits.RData`
  - `data/spec_BLUE_spatial.RData`

- **Agronomic BLUEs:**
  - `data/ag_BLUE_spatial_fits.RData`
  - `data/ag_BLUE_spatial.RData`

- **Malting Quality BLUEss:**
  - `data/mq_BLUP_spatial_fits.RData`
  - `data/mq_BLUP_spatial.RData`

- **Combined VI matrix:**
  - `VImatBLUE.Rdata`

---








### **4. Exploratory Phenotype Analysis**

- **Compute key developmental indicators:**
  - Heading date (HD)
  - Growing degree days (GDD)
  - Julian date (JD) of UAV flights
  - Plot HD, GDD and flights over years 

- **Visualize trait distributions and relationships:**
  - Histograms for AG, MQ, and PCA/FPC traits
  - Correlation matrices and scatter plots

- **Calculating Summary Statistics, CV and overlap**
   -CV for ag, mq and VIs.  

- **Correlation analysis**
   - Ag and MQ x PC scatter with a trend line
   - Ag and MQ x PC heatmap
   - MQ x VI:TP scatter. 
   - MQ c VI:TP heatmap. 
   - VI x PC comp corr




- **Analyze variable importance scores across:**
  - Traits, vegetation indices (VIs), environments, and model types

- Identify informative VI × timepoint combinations
- Select key timepoints for use in multi-trait genomic prediction (MT-GP)


**Outputs:**
- **Summary statistics and correlation tables:**
  - `analysis/CoR_MQ_ag.xlsx`

- **Selected datasets for prediction modeling:**
value
- data/MT_mq_raw
- data/MT_ag_raw
- data/MT_ag_preds







### **5. PCA and FPCA**

- **PCA and FPCA Analysis:**
  - Scree plot, heatmap corr, FPC over time plot.

- **Outputs:**
  - data/F_PCA_scores.Rdata
  - data/F_PCA_pheno_comb.Rdata








### **6. Genomic and Phenomic kernel build. **
- Building relationship matrixes for VIs in different Envs
- Euclidian distance relationship. 

- RAW:
 - ALL different VIs from the specific timepoints 
 - ALL different VIs from all timepoints
 - input: MT_ag_raw and MT_ag_preds
 
 
- BLUE:
 - ALL different VIs from the specific timepoints 
 - ALL different VIs from all timepoints
 - input: VImat and VImatBLUE combined

- **Outputs:**
  -"MT_raw_ wide.Rdata" as MT_raw_w 
  -"MT_pred_wide.Rdata" as MT_pr_w
  - "data/MT_raw_mat.Rdata" as MT_sim_mat_raw
  -  "data/MT_pred_mat.Rdata" as MT_sim_mat_pr
#List of value matrices and relationship matrices for BLUE and raw with combined or selected TP per VIs. 






### **7. Phenomic Prediction mod**



Inputs: 

- data_ag -> pivoted longer
- data_mq -> pivoted longer
- VImat

- ag_BLUE_spatial
- VImatBLUE

- MT_raw_mat
- MT_pred_mat



Outputs: 

LM, RF, LASSO, RR; PCR
- PP_LM_results_ag_raw 
- PP_LM_results_mq_raw 
- PP_LM_results_ag_pred 
- PP_LM_results_mq_pred 
- Var Importance Plots. 



MM with kernel 
- PP_MM_results_ag_pred_comb
- PP_MM_results_ag_pred_sel

- PP_MM_results_ag_raw_comb
- PP_MM_results_ag_raw_sel

- PP_MM_results_MQ_raw_comb
- PP_MM_results_MQ_raw_sel








### **8. Genomic Prediction model build and h2**
- 1 stage
- simple model ang and MQ
- h2 and H2, PAbility

Predict trait BV: 
- GP_MM_results_ag_raw
- GP_MM_results_ag_pred
- GP_MM_results_mq_raw

- GP_PP_results_mq_raw
- GP_PP_results_ag_raw


P_PP_MM_mq
GP_PP_MM_ag

Predict VI pheno, then predict trait BV. 


### **9. MT Genomic Prediction **



### **10. Multi kernel Genomic Prediction **








E1: 
Multi spec pre processing


E2: 
- Run all selected models
- Pick the best performer
- Fit the best model
- Extract predictions
- Set into a dataframe as BLUPS. 


E3. Traits spatial correction model tester
- Predict BLUEs for ag traits with 300 models
- Measure PA for each model prediction in 5 fold CV. 

E44. Trait and  spectral trait correlation analysis
- raw spec cor w MQ traits. 
- 


















# **1. Import data**


**Import agronomic and MQ traits**

```{r}

load("C:/Users/Siim Sepp/NY_WMB/data/ag_BLUE_spatial.RData") # agronomic BLUE
load("C:/Users/Siim Sepp/NY_WMB/data/data_ag.Rdata") # agronomic raw
load("C:/Users/Siim Sepp/NY_WMB/data/mq_BLUP_spatial.RData") # mq BLUP
load("C:/Users/Siim Sepp/NY_WMB/data/data_mq.Rdata") # mq raw
load("C:/Users/Siim Sepp/NY_WMB/data/F_PCA_scores.Rdata") # FPC scres. 
```



**Import raw and predicted VI traits and relationship matrices**


```{r}
load("data/MT_raw_mat.Rdata")
load("data/MT_raw_wide.Rdata")



load("data/MT_pred_mat.Rdata")
load("data/MT_pred_wide.Rdata")
```








- Lot of other options. 
- Does not necessarily need to be a line item of the GPSA 
- Could just do it with a resolution. 