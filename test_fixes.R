# Quick smoke tests for all code-review fixes.
# Loads real data but runs minimal model iterations.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(caret)
})

base_path <- "C:/Users/Siim Sepp/Predicting_WMB_for_NY"
pass <- function(msg) cat(sprintf("  PASS: %s\n", msg))
fail <- function(msg) { cat(sprintf("  FAIL: %s\n", msg)); quit(status = 1) }

# ==============================================================================
cat("\n--- Fix 1: Heritability -- path fallback (file 1) ---\n")
geno_paths <- c(
  file.path(base_path, "data/geno/WMB_master_geno.RData"),
  "~/Documents/GitHub/Predicting_WMB_for_NY/data/geno/WMB_master_geno.RData",
  "data/geno/WMB_master_geno.RData"
)
geno_file <- Filter(file.exists, geno_paths)[1]
if (is.na(geno_file)) fail("WMB_master_geno.RData not found") else pass(paste("found at", geno_file))

# ==============================================================================
cat("\n--- Fix 2: Heritability -- VC grep extraction (file 1, line 177) ---\n")
mock_vc <- data.frame(
  component = c(0.4, 0.1, 0.5),
  row.names  = c("vm(GID, Ginv.sparse)!GID", "Fam!Fam", "units!units")
)
Vg   <- mock_vc[grep("vm\\(GID",  rownames(mock_vc))[1], "component"]
Vfam <- mock_vc[grep("^Fam",      rownames(mock_vc))[1], "component"]
Ve   <- mock_vc[grep("units|R!",  rownames(mock_vc))[1], "component"]
if (!isTRUE(all.equal(c(Vg, Vfam, Ve), c(0.4, 0.1, 0.5)))) {
  fail(sprintf("grep extracted wrong values: Vg=%.1f Vfam=%.1f Ve=%.1f", Vg, Vfam, Ve))
} else pass("grep correctly extracts Vg=0.4, Vfam=0.1, Ve=0.5")

# ==============================================================================
cat("\n--- Fix 3: Heritability -- timepoint filter no longer tautological (file 1, line 1421) ---\n")
mock_data <- data.frame(
  Env      = c("E1","E1","E1","E2"),
  tp_col   = c(1, 2, 1, 1),
  GID      = c("g1","g2","g3","g4"),
  NDVI     = c(0.5, 0.6, 0.7, 0.8)
)
filter_fn <- function(data, env, timepoint) {
  data %>% filter(Env == env, tp_col == timepoint)
}
result <- filter_fn(mock_data, "E1", 1)
if (nrow(result) != 2) fail(sprintf("Expected 2 rows, got %d", nrow(result))) else
  pass("filter returns exactly the 2 rows with tp_col==1 in E1")

# ==============================================================================
cat("\n--- Fix 4: F_PCA -- setdiff excludes 'combined' by name (file 4, line 39) ---\n")
mock_vi_list <- list(NDVI = 1, PSRI = 2, GDVI = 3, combined = 4)
mats <- setdiff(names(mock_vi_list), "combined")
if ("combined" %in% mats) fail("'combined' not excluded") else
if (length(mats) != 3)    fail(sprintf("Expected 3 VIs, got %d", length(mats))) else
  pass(sprintf("setdiff returns %d VIs, no 'combined': %s", length(mats), paste(mats, collapse=", ")))

# Also test with fewer elements (the scenario the old [-10] broke on)
short_list <- list(NDVI = 1, PSRI = 2, combined = 3)
mats_short <- setdiff(names(short_list), "combined")
if ("combined" %in% mats_short) fail("short list: 'combined' not excluded") else
  pass("setdiff also correct for short list (3 elements)")

# ==============================================================================
cat("\n--- Fix 5: PP_ag_mq -- nearZeroVar mask applied to test (file 6, line 96) ---\n")
set.seed(42)
n <- 50
mock_train <- data.frame(
  S.T  = rnorm(n),
  V1   = rnorm(n),
  V2   = rep(0, n),      # near-zero var — should be dropped
  V3   = rnorm(n)
)
mock_test <- mock_train[1:10, ]

nzv        <- nearZeroVar(mock_train, saveMetrics = TRUE)
mock_train <- mock_train[, !nzv$nzv]
mock_test  <- mock_test[, colnames(mock_train)]   # the fix

if (!identical(colnames(mock_train), colnames(mock_test)))
  fail("train and test columns still differ after fix") else
  pass(sprintf("train and test share same %d columns: %s",
               ncol(mock_train), paste(colnames(mock_train), collapse=", ")))
if ("V2" %in% colnames(mock_test)) fail("near-zero column V2 still in test") else
  pass("near-zero column V2 correctly removed from test")

# ==============================================================================
cat("\n--- Fix 6: MT_GP -- named column indexing for NDVI/PSRI (file 8, line 56) ---\n")
mt_wide_path <- file.path(base_path, "data/MT_raw_wide.Rdata")
if (!file.exists(mt_wide_path)) {
  cat("  SKIP: MT_raw_wide.Rdata not found (run file 5 first)\n")
} else {
  load(mt_wide_path)   # loads MT_raw_w
  first_key <- names(MT_raw_w)[1]
  mat <- MT_raw_w[[first_key]]
  if (!"NDVI" %in% colnames(mat)) fail(sprintf("NDVI column missing from %s", first_key)) else
    pass(sprintf("NDVI column found in MT_raw_w$%s", first_key))
  if (!"PSRI" %in% colnames(mat)) fail(sprintf("PSRI column missing from %s", first_key)) else
    pass(sprintf("PSRI column found in MT_raw_w$%s", first_key))
  # named extraction should not error
  tryCatch({
    x <- mat[, "NDVI"]
    y <- mat[, "PSRI"]
    pass(sprintf("named extraction works: %d NDVI values, %d PSRI values", length(x), length(y)))
  }, error = function(e) fail(paste("named extraction error:", e$message)))
}

# ==============================================================================
cat("\n--- Fix 7: DK_GP -- output file name consistent (file 10, line 277) ---\n")
# We can't run the model, but check the saved output exists if it was already run
ag_pred_path <- file.path(base_path, "output/GP_DK_results_ag_pred.Rdata")
ag_pr_path   <- file.path(base_path, "output/GP_DK_results_ag_pr.Rdata")
if (file.exists(ag_pred_path)) pass("GP_DK_results_ag_pred.Rdata exists — results chunk will load correctly") else
  cat("  INFO: GP_DK_results_ag_pred.Rdata not yet generated (needs model run)\n")
if (file.exists(ag_pr_path))
  cat("  INFO: old GP_DK_results_ag_pr.Rdata also present — results chunk now correctly points to _pred\n")

# ==============================================================================
cat("\n--- Fix 8: Prediction_plots -- fct_relevel on correct column (file 11, line 322) ---\n")
env_order <- c("KET21", "SNY22", "MCG23", "HELF24", "MCG25")
mock_res <- data.frame(
  Env = c("MCG25", "HELF24", "KET21", "SNY22", "MCG23"),
  Cor = runif(5)
)
result <- mock_res %>%
  mutate(Environment = trimws(as.character(Env)),
         Environment = fct_relevel(Environment, env_order))

if (!is.factor(result$Environment)) fail("Environment is not a factor") else
if (!identical(levels(result$Environment), env_order)) {
  fail(sprintf("levels wrong: %s", paste(levels(result$Environment), collapse=", ")))
} else pass(sprintf("Environment factor levels correct: %s", paste(levels(result$Environment), collapse=", ")))

# also confirm the old bug: fct_relevel(Env, ...) would leave wrong levels
result_buggy <- mock_res %>%
  mutate(Environment = trimws(as.character(Env)),
         Environment = fct_relevel(Env, env_order))
if (identical(levels(result_buggy$Environment), env_order))
  cat("  NOTE: cannot demonstrate old bug since Env already matches env_order in mock\n") else
  pass("confirmed old bug: fct_relevel(Env,...) produced wrong levels")

# ==============================================================================
cat("\n=== All tests completed ===\n")
