# Pipeline smoke test.
# Purls each Rmd, reduces iteration counts to minimum, sources in sequence.
# Each file's saved outputs land in the project directories as normal,
# so downstream files can load them.

suppressPackageStartupMessages(library(knitr))
Sys.setenv(RSTUDIO_PANDOC =
  "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

base_path <- "C:/Users/Siim Sepp/Predicting_WMB_for_NY"
setwd(base_path)

results <- list()

run_rmd <- function(label, filename, patches = character(0)) {
  cat(sprintf("\n==============================\n[%s] %s\n", label, filename))
  rmd_path <- file.path(base_path, "analysis", filename)

  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp), add = TRUE)

  knitr::purl(rmd_path, output = tmp, quiet = TRUE, documentation = 0)
  code <- readLines(tmp)

  for (i in seq_along(patches)) {
    code <- gsub(names(patches)[i], patches[i], code, perl = TRUE)
  }

  writeLines(code, tmp)

  tryCatch({
    source(tmp, local = FALSE, echo = FALSE)   # local=FALSE -> globalenv
    cat("  --> PASS\n")
    results[[label]] <<- "PASS"
  }, error = function(e) {
    cat(sprintf("  --> FAIL: %s\n", e$message))
    results[[label]] <<- paste("FAIL:", e$message)
  })
}

# ── shared speed patches ───────────────────────────────────────────────────────
speed <- c(
  "(?<=(reps\\s{0,5}=\\s{0,5}))\\d+"   = "1",
  "(?<=(folds\\s{0,5}=\\s{0,5}))\\d+"  = "2",
  "(?<=(cycles\\s{0,5}<-\\s{0,5}))\\d+"= "1",
  "(?<=(reps\\s{0,5}<-\\s{0,5}))\\d+"  = "1"
)

# ── per-file environment limits ────────────────────────────────────────────────
# Restrict envs/traits to first element so model loops finish in seconds.
env1  <- c("(?<=envs\\s{0,5}<-\\s{0,5})unique.*"              = 'unique(data_ag$Env)[1]',
           "(?<=envs\\s{0,5}<-\\s{0,5})levels.*"             = 'levels(data_mq$Env)[1]',
           "(?<=envs\\s{0,5}=\\s{0,5})unique.*"              = 'unique(data$Env)[1]')
trait1 <- c("(?<=traits\\s{0,5}<-\\s{0,5})colnames.*"        = 'colnames(data_ag)[7]')

# ── run each file in order ─────────────────────────────────────────────────────

# Files 1–5 do single REML fits or data prep — run at full scale so their
# saved outputs (ag_BLUE_spatial.RData, MT_pred_mat.Rdata, etc.) are complete
# and usable by the downstream CV files.
run_rmd("01_heritability", "1. Trait heritability analysis.Rmd")

run_rmd("02_exploratory",  "2. Exploratory_pheno_ analysis.Rmd")

run_rmd("03_spectral_BLUE", "3. Spectral_data_BLUE.Rmd")

run_rmd("04_FPCA",  "4. F_PCA.Rmd")

run_rmd("05_kernels", "5. Phenomic kernel_build.Rmd")

run_rmd("06_PP",    "6. PP_ag_mq.Rmd",   speed)

run_rmd("07_GP",    "7. GP_ag_mq.Rmd",   speed)

run_rmd("08_MT_GP", "8. MT_GP.Rmd",      speed)

run_rmd("10_DK_GP", "10.DK_GP.Rmd",      speed)

run_rmd("11_plots", "11. Prediction_plots.Rmd")

# ── summary ───────────────────────────────────────────────────────────────────
cat("\n\n============================== SUMMARY ==============================\n")
for (lbl in names(results)) {
  status <- results[[lbl]]
  cat(sprintf("  %-20s  %s\n", lbl, status))
}
passes <- sum(sapply(results, function(x) x == "PASS"))
cat(sprintf("\n%d / %d files passed.\n", passes, length(results)))
