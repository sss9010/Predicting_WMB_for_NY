library(readxl)

# ---- Build mapping ----
map <- read_excel("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/map_UMN.xlsx")

# Normalize: strip leading zeros from number part (BAROREI0011 -> BAROREI11)
normalize_id <- function(ids) sub("^(BAROREI)0*(\\d+)$", "\\1\\2", ids)

map$Entry_norm <- normalize_id(map$Entry)
lookup <- setNames(map$GID, map$Entry_norm)  # named vector: BAROREI11 -> 25BS-NAM-xxx

cat("Mapping built:", nrow(map), "entries\n")
cat("Example: BAROREI11 ->", lookup["BAROREI11"], "\n\n")

# ---- Helper: rename a vector of sample IDs ----
rename_ids <- function(ids) {
  norms <- normalize_id(ids)
  ifelse(norms %in% names(lookup), lookup[norms], ids)
}

# ---- Helper: rename IDs with .suffix (FDT/mFDT) ----
rename_dot_ids <- function(cols) {
  sapply(cols, function(col) {
    m <- regmatches(col, regexpr("^(BAROREI\\d+)(\\..*)?$", col))
    if (length(m) == 0) return(col)
    base <- sub("\\..*$", "", col)
    suffix <- sub("^BAROREI\\d+", "", col)
    norm_base <- normalize_id(base)
    new_base <- if (norm_base %in% names(lookup)) lookup[norm_base] else base
    paste0(new_base, suffix)
  }, USE.NAMES = FALSE)
}

# ---- Helper: rename id{n}.BAROREI{n} columns (IDnum VCF) ----
rename_idnum_ids <- function(cols) {
  sapply(cols, function(col) {
    m <- regmatches(col, regexpr("^(id\\d+\\.)(BAROREI\\d+)$", col))
    if (length(m) == 0) return(col)
    prefix <- sub("(id\\d+\\.)(BAROREI\\d+)$", "\\1", col)
    barorei <- sub("^id\\d+\\.", "", col)
    norm_b <- normalize_id(barorei)
    new_b <- if (norm_b %in% names(lookup)) lookup[norm_b] else barorei
    paste0(prefix, new_b)
  }, USE.NAMES = FALSE)
}

# ---- Helper: rewrite a large file - only modify the header line ----
rewrite_header_line <- function(filepath, modify_fn, is_vcf = FALSE, n_fixed_cols = 0) {
  tmpfile <- paste0(filepath, ".tmp")

  con_in  <- file(filepath, open = "r", encoding = "UTF-8")
  con_out <- file(tmpfile,  open = "w", encoding = "UTF-8")

  modified <- FALSE
  chunk_size <- 50000

  tryCatch({
    repeat {
      lines <- readLines(con_in, n = chunk_size, warn = FALSE)
      if (length(lines) == 0) break

      if (!modified) {
        if (is_vcf) {
          # Pass through ## comment lines unchanged, modify #CHROM line
          for (i in seq_along(lines)) {
            if (grepl("^#CHROM", lines[i])) {
              cols <- strsplit(lines[i], "\t")[[1]]
              fixed <- cols[1:9]
              sample_cols <- cols[10:length(cols)]
              new_samples <- modify_fn(sample_cols)
              lines[i] <- paste(c(fixed, new_samples), collapse = "\t")
              modified <- TRUE
              break
            }
          }
          # If #CHROM not found yet in this chunk, pass through and continue
          writeLines(lines, con_out)
        } else {
          # First line is header
          cols <- strsplit(lines[1], "\t")[[1]]
          fixed <- if (n_fixed_cols > 0) cols[1:n_fixed_cols] else character(0)
          sample_cols <- cols[(n_fixed_cols + 1):length(cols)]
          new_samples <- modify_fn(sample_cols)
          lines[1] <- paste(c(fixed, new_samples), collapse = "\t")
          modified <- TRUE
          writeLines(lines, con_out)
        }
      } else {
        writeLines(lines, con_out)
      }
    }
  }, finally = {
    close(con_in)
    close(con_out)
  })

  if (!modified) {
    file.remove(tmpfile)
    warning("No header found to modify in: ", filepath)
    return(invisible(FALSE))
  }

  file.remove(filepath)
  file.rename(tmpfile, filepath)
  invisible(TRUE)
}

# ---- Helper: rewrite Samples.txt (rename Sample.ID column values) ----
rewrite_samples_txt <- function(filepath) {
  df <- read.delim(filepath, check.names = FALSE, stringsAsFactors = FALSE)
  # Rename Sample.ID column if present
  id_col <- intersect(c("Sample.ID", "Sample ID"), names(df))
  if (length(id_col) == 0) {
    cat("  WARNING: no Sample.ID column found in", basename(filepath), "\n")
    return(invisible(FALSE))
  }
  df[[id_col]] <- rename_ids(df[[id_col]])
  write.table(df, filepath, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(TRUE)
}

# ---- Helper: rewrite VCF-Samples.txt (tab-delimited, rename ID columns) ----
rewrite_vcf_samples_txt <- function(filepath) {
  if (file.info(filepath)$size == 0) { cat("  (empty — skipped)\n"); return(invisible(FALSE)) }
  df <- read.delim(filepath, check.names = FALSE, stringsAsFactors = FALSE)
  cat("  Columns:", paste(names(df), collapse=", "), "\n")

  # Rename Sample.ID / Sample ID column
  for (col in c("Sample.ID", "Sample ID")) {
    if (col %in% names(df)) df[[col]] <- rename_ids(df[[col]])
  }
  # Rename VCF Name column (plain BAROREI)
  for (col in c("VCF.Name", "VCF Name")) {
    if (col %in% names(df)) df[[col]] <- rename_ids(df[[col]])
  }
  # Rename IDnumVCF Name column (id{n}.BAROREI{n} format)
  for (col in c("IDnumVCF.Name", "IDnumVCF Name")) {
    if (col %in% names(df)) df[[col]] <- rename_idnum_ids(df[[col]])
  }
  write.table(df, filepath, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(TRUE)
}

# ============================================================
# Process all 3K files
# ============================================================
base <- "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K"

# ---- 1. VCF files (.MorexV3.vcf) ----
vcf_files <- c(
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri.MorexV3.vcf"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO.MorexV3.vcf")
)
for (f in vcf_files) {
  cat("Processing VCF:", basename(f), "...")
  ok <- rewrite_header_line(f, rename_ids, is_vcf = TRUE)
  cat(if (ok) "DONE\n" else "SKIPPED\n")
}

# ---- 2. IDnum VCF files (.MorexV3.IDnum.vcf) ----
idnum_files <- c(
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri.MorexV3.IDnum.vcf"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO.MorexV3.IDnum.vcf")
)
for (f in idnum_files) {
  cat("Processing IDnum VCF:", basename(f), "...")
  ok <- rewrite_header_line(f, rename_idnum_ids, is_vcf = TRUE)
  cat(if (ok) "DONE\n" else "SKIPPED\n")
}

# ---- 3. AB files (header: 3 fixed cols: Name, Chr, Pos) ----
ab_files <- c(
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-AB.txt"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-AB.txt"),
  file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-AB.txt")
)
for (f in ab_files) {
  cat("Processing AB:", basename(f), "...")
  ok <- rewrite_header_line(f, rename_ids, n_fixed_cols = 3)
  cat(if (ok) "DONE\n" else "SKIPPED\n")
}

# ---- 4. FDT files (header: 5 fixed cols, then BAROREI.suffix) ----
fdt_files <- c(
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-FDT.txt"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-FDT.txt"),
  file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-FDT.txt")
)
for (f in fdt_files) {
  if (!file.exists(f)) { cat("FDT not found:", basename(f), "\n"); next }
  cat("Processing FDT:", basename(f), "...")
  ok <- rewrite_header_line(f, rename_dot_ids, n_fixed_cols = 5)
  cat(if (ok) "DONE\n" else "SKIPPED\n")
}

# ---- 5. mFDT files (header: 5 fixed cols, then BAROREI.suffix) ----
mfdt_files <- c(
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-mFDT.txt"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-mFDT.txt")
)
for (f in mfdt_files) {
  if (!file.exists(f)) { cat("mFDT not found:", basename(f), "\n"); next }
  cat("Processing mFDT:", basename(f), "...")
  ok <- rewrite_header_line(f, rename_dot_ids, n_fixed_cols = 5)
  cat(if (ok) "DONE\n" else "SKIPPED\n")
}

# ---- 6. R files (header: 2 fixed cols: CHROM, POS) ----
r_files <- c(
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-R.txt"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-R.txt"),
  file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-R.txt")
)
for (f in r_files) {
  if (!file.exists(f)) { cat("R file not found:", basename(f), "\n"); next }
  cat("Processing R:", basename(f), "...")
  ok <- rewrite_header_line(f, rename_ids, n_fixed_cols = 2)
  cat(if (ok) "DONE\n" else "SKIPPED\n")
}

# ---- 7. Samples.txt files ----
samples_files <- c(
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-Samples.txt"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-Sample.txt"),
  file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-Samples.txt")
)
for (f in samples_files) {
  if (!file.exists(f)) { cat("Samples file not found:", basename(f), "\n"); next }
  cat("Processing Samples:", basename(f), "...")
  ok <- rewrite_samples_txt(f)
  cat(if (ok) "DONE\n" else "SKIPPED\n")
}

# ---- 8. VCF-Samples and poorQC-Samples files ----
vcf_sample_files <- c(
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-VCF-Samples.txt"),
  file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-poorQC_VCF-Samples.txt"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-Sample.txtVCF-Samples.txt"),
  file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-poorQC_VCF-Samples.txt"),
  file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-VCF-Samples.txt"),
  file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-poorQC_VCF-Samples.txt")
)
for (f in vcf_sample_files) {
  if (!file.exists(f)) { cat("VCF-Samples file not found:", basename(f), "\n"); next }
  cat("Processing VCF-Samples:", basename(f), "\n")
  ok <- rewrite_vcf_samples_txt(f)
  cat(if (ok) "  DONE\n" else "  SKIPPED\n")
}

cat("\n=== All files processed ===\n")

# ---- Verification: spot-check a few files ----
cat("\n--- Verification spot-check ---\n")

# Check VCF header
con <- file(file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri.MorexV3.vcf"), "r")
lines <- readLines(con, n=20)
close(con)
chrom_line <- lines[grepl("^#CHROM", lines)]
cols <- strsplit(chrom_line, "\t")[[1]]
cat("VCF Tri - first 15 sample cols:\n")
print(cols[10:min(24, length(cols))])

# Check AB header
con2 <- file(file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-AB.txt"), "r")
h2 <- readLines(con2, n=1)
close(con2)
cols2 <- strsplit(h2, "\t")[[1]]
cat("\nAB COMBO - first 15 sample cols:\n")
print(cols2[4:min(18, length(cols2))])
