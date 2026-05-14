base <- "C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_X"

# VCF header
con <- file(file.path(base, "UWMad2024Gutierrez-B3K_01-06.MorexV3.vcf"), "r")
lines <- readLines(con, n=30); close(con)
chrom <- lines[grepl("^#CHROM", lines)]
cols <- strsplit(chrom, "\t")[[1]]
cat("VCF - N samples:", length(cols)-9, "\n")
cat("First 8 sample IDs:", paste(cols[10:17], collapse=", "), "\n\n")

# AB header
con2 <- file(file.path(base, "UWMad2024Gutierrez-B3K_01-06-AB.txt"), "r")
h <- readLines(con2, n=1); close(con2)
cols2 <- strsplit(h, "\t")[[1]]
cat("AB - N samples:", length(cols2)-3, "\n")
cat("First 8 sample IDs:", paste(cols2[4:11], collapse=", "), "\n\n")

# CSV
csv <- read.csv(file.path(base, "UWMad2024Gutierrez-B3K_01-06.csv"), nrows=5)
cat("CSV columns:", paste(names(csv), collapse=", "), "\n")
cat("First 5 rows:\n")
print(head(csv, 5))

# Samples
s <- read.delim(file.path(base, "UWMad2024Gutierrez-B3K_01-06-Samples.txt"))
cat("\nSamples columns:", paste(names(s), collapse=", "), "\n")
cat("First 5 Sample IDs:", paste(head(s$Sample.ID, 5), collapse=", "), "\n")
cat("Total samples:", nrow(s), "\n")
