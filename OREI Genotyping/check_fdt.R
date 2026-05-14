con <- file("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-FDT.txt", "r")
header <- readLines(con, n=1)
close(con)
cols <- strsplit(header, "\t")[[1]]
cat("Total columns:", length(cols), "\n")
cat("First 15 columns:\n")
print(cols[1:15])
cat("...\nLast 5 columns:\n")
print(tail(cols, 5))

# Check mFDT
con2 <- file("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-mFDT.txt", "r")
h2 <- readLines(con2, n=1)
close(con2)
cols2 <- strsplit(h2, "\t")[[1]]
cat("\nmFDT total cols:", length(cols2), "\n")
cat("First 10:\n")
print(cols2[1:10])

# Check R file
con3 <- file("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-R.txt", "r")
h3 <- readLines(con3, n=1)
close(con3)
cols3 <- strsplit(h3, "\t")[[1]]
cat("\nR file total cols:", length(cols3), "\n")
cat("First 10:\n")
print(cols3[1:10])

# Check Samples.txt
s <- read.delim("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-Samples.txt")
cat("\nSamples.txt columns:", paste(names(s), collapse=", "), "\n")
cat("First 5 rows:\n")
print(head(s[,1:5], 5))

# Check VCF-Samples.txt
vk <- readLines("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri-VCF-Samples.txt", n=10)
cat("\nVCF-Samples (first 10):\n")
print(vk)
