library(readxl)

map <- read_excel("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/map_UMN.xlsx")
cat("Map BAROREI IDs (first 10):", paste(head(map$Entry, 10), collapse=", "), "\n")
cat("Map GID IDs (first 10):", paste(head(map$GID, 10), collapse=", "), "\n\n")

# Check VCF headers
vcf_files <- c(
  "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri.MorexV3.vcf",
  "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO.MorexV3.vcf"
)

for (f in vcf_files) {
  lines <- readLines(f, n=50)
  header_line <- lines[grepl("^#CHROM", lines)]
  cols <- strsplit(header_line, "\t")[[1]]
  sample_cols <- cols[10:length(cols)]
  cat("File:", basename(f), "\n")
  cat("  N samples:", length(sample_cols), "\n")
  cat("  First 10:", paste(head(sample_cols, 10), collapse=", "), "\n")
  cat("  Last 5:", paste(tail(sample_cols, 5), collapse=", "), "\n")
  in_map <- sample_cols[sample_cols %in% map$Entry]
  cat("  Matching map entries:", length(in_map), "\n\n")
}

# Also check AB files
ab_files <- c(
  "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-AB.txt",
  "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-AB.txt",
  "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-14-COMBO/OREI2024-B3K_01-14-AB.txt"
)

for (f in ab_files) {
  if (!file.exists(f)) { cat("NOT FOUND:", f, "\n"); next }
  header <- strsplit(readLines(f, n=1), "\t")[[1]]
  sample_cols <- header[4:length(header)]
  cat("AB File:", basename(f), "\n")
  cat("  N samples:", length(sample_cols), "\n")
  cat("  First 10:", paste(head(sample_cols, 10), collapse=", "), "\n")
  cat("  Last 5:", paste(tail(sample_cols, 5), collapse=", "), "\n")
  in_map <- sample_cols[sample_cols %in% map$Entry]
  cat("  Matching map entries:", length(in_map), "\n\n")
}
