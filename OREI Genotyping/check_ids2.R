library(readxl)

map <- read_excel("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/map_UMN.xlsx")
cat("All map Entry IDs:\n")
print(sort(map$Entry))
cat("\n")

# Check which DualO samples match
vcf <- "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO.MorexV3.vcf"
lines <- readLines(vcf, n=50)
header_line <- lines[grepl("^#CHROM", lines)]
cols <- strsplit(header_line, "\t")[[1]]
sample_cols <- cols[10:length(cols)]
in_map <- sample_cols[sample_cols %in% map$Entry]
cat("DualO samples in map:\n")
print(sort(in_map))
cat("\n")

# Check samples NOT in map (from DualO)
not_in_map <- sample_cols[!sample_cols %in% map$Entry]
cat("DualO samples NOT in map (first 20):\n")
print(head(sort(not_in_map), 20))

# Also check the Samples.txt file for DualO to see the actual sample IDs
samples_file <- "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_11-14_DualO/OREI2024-B3K_11-14_DualO-Sample.txt"
if (file.exists(samples_file)) {
  s <- read.delim(samples_file)
  cat("\nSamples file columns:", paste(names(s), collapse=", "), "\n")
  print(head(s, 10))
}
