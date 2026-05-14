con <- file("C:/Users/Siim Sepp/NY_WMBCL/OREI Genotyping/OREI_Genotyping_3K/OREI2024-B3K_01-10_Tri/OREI2024-B3K_01-10_Tri.MorexV3.IDnum.vcf", "r")
lines <- readLines(con, n=30)
close(con)
# Find #CHROM line
chrom_line <- lines[grepl("^#CHROM", lines)]
cols <- strsplit(chrom_line, "\t")[[1]]
cat("IDnum VCF - first 15 columns:\n")
print(cols[1:15])
cat("Last 3:\n")
print(tail(cols, 3))
