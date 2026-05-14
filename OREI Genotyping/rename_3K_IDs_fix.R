library(readxl)

map <- read_excel("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/map_UMN.xlsx")
normalize_id <- function(ids) sub("^(BAROREI)0*(\\d+)$", "\\1\\2", ids)
map$Entry_norm <- normalize_id(map$Entry)
lookup <- setNames(map$GID, map$Entry_norm)

rename_ids <- function(ids) {
  norms <- normalize_id(ids)
  ifelse(norms %in% names(lookup), lookup[norms], ids)
}
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

rewrite_vcf_samples_txt <- function(filepath) {
  if (file.info(filepath)$size == 0) {
    cat("  Skipping empty file\n")
    return(invisible(FALSE))
  }
  df <- read.delim(filepath, check.names = FALSE, stringsAsFactors = FALSE)
  cat("  Columns:", paste(names(df), collapse=", "), "\n")
  for (col in c("Sample.ID", "Sample ID")) {
    if (col %in% names(df)) df[[col]] <- rename_ids(df[[col]])
  }
  for (col in c("VCF.Name", "VCF Name")) {
    if (col %in% names(df)) df[[col]] <- rename_ids(df[[col]])
  }
  for (col in c("IDnumVCF.Name", "IDnumVCF Name")) {
    if (col %in% names(df)) df[[col]] <- rename_idnum_ids(df[[col]])
  }
  write.table(df, filepath, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(TRUE)
}

base <- "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K"

remaining <- c(
  file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-VCF-Samples.txt"),
  file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-poorQC_VCF-Samples.txt")
)

for (f in remaining) {
  cat("Processing:", basename(f), "\n")
  ok <- rewrite_vcf_samples_txt(f)
  cat(if (isTRUE(ok)) "  DONE\n" else "  SKIPPED\n")
}

cat("\n=== Verification: spot-check renamed IDs ===\n")

# Check VCF Tri header
con <- file(file.path(base, "OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri.MorexV3.vcf"), "r")
lines <- readLines(con, n=20); close(con)
chrom_line <- lines[grepl("^#CHROM", lines)]
cols <- strsplit(chrom_line, "\t")[[1]]
cat("VCF Tri - sample cols 10-24:\n")
print(cols[10:min(24, length(cols))])

# Check DualO VCF
con2 <- file(file.path(base, "OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO.MorexV3.vcf"), "r")
lines2 <- readLines(con2, n=20); close(con2)
chrom2 <- lines2[grepl("^#CHROM", lines2)]
cols2 <- strsplit(chrom2, "\t")[[1]]
cat("\nVCF DualO - first 15 sample cols:\n")
print(cols2[10:min(24, length(cols2))])

# Check AB COMBO
con3 <- file(file.path(base, "OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-AB.txt"), "r")
h3 <- readLines(con3, n=1); close(con3)
cols3 <- strsplit(h3, "\t")[[1]]
cat("\nAB COMBO - sample cols 4-18:\n")
print(cols3[4:min(18, length(cols3))])
