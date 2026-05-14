library(readxl)

map <- read_excel("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/map_UMN.xlsx")

# Normalize: strip leading zeros from number part (BAROREI0011 -> BAROREI11)
normalize_id <- function(ids) {
  sub("^(BAROREI)0*(\\d+)$", "\\1\\2", ids)
}

map$Entry_norm <- normalize_id(map$Entry)
cat("Normalized map entries (first 10):", paste(head(map$Entry_norm, 10), collapse=", "), "\n\n")

# Check all 3K files
files_to_check <- list(
  vcf_tri    = "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri.MorexV3.vcf",
  vcf_dualo  = "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO.MorexV3.vcf",
  ab_tri     = "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-AB.txt",
  ab_dualo   = "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-AB.txt",
  ab_combo   = "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-AB.txt"
)

get_header_samples <- function(f, skip_cols=9) {
  line1 <- readLines(f, n=1)
  # Skip comment lines for VCF
  if (grepl("^##", line1)) {
    lines <- readLines(f, n=200)
    line1 <- lines[grepl("^#CHROM", lines)]
  }
  cols <- strsplit(line1, "\t")[[1]]
  cols[(skip_cols+1):length(cols)]
}

for (nm in names(files_to_check)) {
  f <- files_to_check[[nm]]
  if (!file.exists(f)) { cat(nm, ": FILE NOT FOUND\n"); next }
  skip <- if (grepl("vcf$", f)) 9 else 3
  samples <- get_header_samples(f, skip)
  samples_norm <- normalize_id(samples)
  in_map <- samples_norm[samples_norm %in% map$Entry_norm]
  cat(nm, "(", length(samples), "samples ):", length(in_map), "match map\n")
}
