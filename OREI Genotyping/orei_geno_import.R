# OREI 2024 B3K genotyping – import and QC preparation
# Uses the COMBO batch (plates 01-14) as the primary dataset.
# Outputs: marker_map, geno_mat (numeric 0/1/2/NA), marker_qc, sample_qc

library(data.table)

# ── 1. Paths ──────────────────────────────────────────────────────────────────

batch_dir <- file.path("OREI Genotyping", "OREI2024-B3K_01-14-COMBO")
prefix    <- "OREI2024-B3K_01-14"

path <- function(suffix) file.path(batch_dir, paste0(prefix, suffix))

# ── 2. Sample metadata ────────────────────────────────────────────────────────

samples_meta <- fread(path("-Samples.txt"),         sep = "\t", check.names = FALSE)
pass_meta    <- fread(path("-VCF-Samples.txt"),     sep = "\t", check.names = FALSE)
fail_meta    <- fread(path("-poorQC_VCF-Samples.txt"), sep = "\t", check.names = FALSE)

# Remove blank trailing rows that GenomeStudio sometimes writes
fail_meta <- fail_meta[!is.na(`Sample ID`) & nzchar(`Sample ID`)]

samples_meta[, qc_pass := `Sample ID` %in% pass_meta[["Sample ID"]]
             & !(`Sample ID` %in% fail_meta[["Sample ID"]])]

# ── 3. AB genotype calls ──────────────────────────────────────────────────────
# comment.char="" prevents fread treating the leading "#" as a comment marker

ab <- fread(path("-AB.txt"), sep = "\t", header = TRUE)#,
            #comment.char = "", check.names = FALSE)

setnames(ab, names(ab)[1], "marker_id")   # "#ID" → "marker_id"

marker_map <- ab[, .(marker_id,
                     chrom = CHROM,
                     pos   = as.integer(POS))]

sample_ids <- setdiff(names(ab), c("marker_id", "CHROM", "POS"))

# ── 4. Encode AA=0, AB=1, BB=2, NC=NA ────────────────────────────────────────

ab_to_int <- function(m) {
  out        <- matrix(NA_integer_, nrow(m), ncol(m), dimnames = dimnames(m))
  out[m == "AA"] <- 0L
  out[m == "AB"] <- 1L
  out[m == "BB"] <- 2L
  out
}

char_mat <- as.matrix(ab[, ..sample_ids])
rownames(char_mat) <- ab$marker_id

geno_mat <- ab_to_int(char_mat)   # markers × samples, integer 0/1/2/NA

# ── 5. Marker-level QC ───────────────────────────────────────────────────────

n_samples   <- ncol(geno_mat)
miss_marker <- rowMeans(is.na(geno_mat))

allele_freq <- rowMeans(geno_mat, na.rm = TRUE) / 2   # frequency of B allele
maf         <- pmin(allele_freq, 1 - allele_freq)

obs_het_marker <- rowMeans(geno_mat == 1L, na.rm = TRUE)

marker_qc <- data.table(
  marker_id = marker_map$marker_id,
  chrom     = marker_map$chrom,
  pos       = marker_map$pos,
  miss_rate = miss_marker,
  maf       = maf,
  obs_het   = obs_het_marker
)

# ── 6. Sample-level QC ───────────────────────────────────────────────────────

miss_sample <- colMeans(is.na(geno_mat))
obs_het_sample <- colMeans(geno_mat == 1L, na.rm = TRUE)

meta_idx <- match(sample_ids, samples_meta[["Sample ID"]])

sample_qc <- data.table(
  sample_id       = sample_ids,
  plate           = samples_meta[["Sample Plate"]][meta_idx],
  call_rate_array = samples_meta[["Call Rate"]][meta_idx],
  qc_pass_array   = samples_meta[["qc_pass"]][meta_idx],
  call_rate_geno  = 1 - miss_sample,
  obs_het         = obs_het_sample
)

# ── 7. Summary ────────────────────────────────────────────────────────────────

cat(sprintf("Markers  : %d\n",                             nrow(marker_qc)))
cat(sprintf("Samples  : %d\n",                             nrow(sample_qc)))
cat(sprintf("  Pass array QC    : %d\n",                   sum(sample_qc$qc_pass_array, na.rm = TRUE)))
cat(sprintf("  Fail array QC    : %d\n",                   sum(!sample_qc$qc_pass_array, na.rm = TRUE)))
cat(sprintf("Markers >20%% miss  : %d\n",                  sum(marker_qc$miss_rate > 0.20)))
cat(sprintf("Markers MAF <0.05  : %d\n",                   sum(marker_qc$maf < 0.05, na.rm = TRUE)))
cat(sprintf("Samples call <0.80 : %d\n",                   sum(sample_qc$call_rate_geno < 0.80)))

# ── Objects available for downstream QC ──────────────────────────────────────
#   geno_mat    – integer matrix, markers × samples (0/1/2/NA)
#   marker_map  – data.table: marker_id / chrom / pos
#   marker_qc   – per-marker: miss_rate, maf, obs_het
#   sample_qc   – per-sample: call rates (array + recalculated), obs_het, QC flag
#   samples_meta – full GenomeStudio sample metadata
