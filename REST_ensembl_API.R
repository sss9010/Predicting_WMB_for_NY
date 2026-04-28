
# ── Libraries ─────────────────────────────────────────────────────────────────
library(httr)
library(jsonlite)
library(biomaRt)
library(dplyr)
library(readr)

# ══════════════════════════════════════════════════════════════════════════════
# HELPER: parse JSON gene response from Ensembl REST
# ══════════════════════════════════════════════════════════════════════════════
parse_ensembl_genes <- function(response) {
  genes <- fromJSON(content(response, "text"))
  if (length(genes) == 0) return(data.frame(
    gene_id = character(), gene_start = integer(), gene_end = integer(),
    strand = integer(), biotype = character(), canonical_transcript = character(),
    stringsAsFactors = FALSE))
  data.frame(
    gene_id              = genes$gene_id,
    gene_start           = genes$start,
    gene_end             = genes$end,
    strand               = genes$strand,
    biotype              = genes$biotype,
    canonical_transcript = genes$canonical_transcript,
    stringsAsFactors = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION: pull candidate genes within a flanking window of each SNP
# Arguments:
#   SNPs        - data frame with chromosome (col 1) and position (col 2)
#   species     - Ensembl species name (Hordeum_vulgare / Triticum_aestivum)
#   flank_length- bp upstream and downstream of each SNP to query
# Returns: candidate_genes data frame
# ══════════════════════════════════════════════════════════════════════════════
extract_candidates_snp <- function(SNPs, species, flank_length) {
  candidate_genes <- data.frame(
    neighboring_snp = integer(), gene_id = character(),
    gene_start = integer(), gene_end = integer(),
    strand = integer(), biotype = character(),
    canonical_transcript = character(),
    stringsAsFactors = FALSE)

  server    <- "https://rest.ensembl.org"
  CHUNK_MAX <- 4999999L   # Ensembl REST inclusive-coord limit: end-start+1 <= 5000000
  # Ensembl Plants barley chromosomes are named "1H", "2H", etc.
  # SNP Chr column contains the number only (e.g. 6); append "H" for barley.
  chrom_suffix <- if (grepl("hordeum_vulgare", tolower(species))) "H" else ""

  for (i in seq_len(nrow(SNPs))) {
    current_chrom <- paste0(SNPs[i, 1], chrom_suffix)
    current_snp   <- as.integer(SNPs[i, 2])
    start_bp <- max(1L, current_snp - as.integer(flank_length))
    end_bp   <- current_snp + as.integer(flank_length)

    # Split into ≤5 Mb chunks if the flanking window is too large
    chunk_starts <- seq(start_bp, end_bp, by = CHUNK_MAX)
    chunk_ends   <- pmin(chunk_starts + CHUNK_MAX - 1L, end_bp)

    for (j in seq_along(chunk_starts)) {
      ext <- paste0("/overlap/region/", tolower(species), "/",
                    current_chrom, ":",
                    chunk_starts[j], "-", chunk_ends[j],
                    "?feature=gene")
      message("  Querying: ", server, ext)
      request <- GET(paste0(server, ext), accept("application/json"))
      stop_for_status(request)

      parsed_genes <- parse_ensembl_genes(request)
      if (nrow(parsed_genes) > 0)
        candidate_genes <- rbind(candidate_genes,
                                 cbind(neighboring_snp = current_snp, parsed_genes))
    }
  }
  # Remove duplicates from genes spanning chunk boundaries
  candidate_genes <- dplyr::distinct(candidate_genes)
  invisible(candidate_genes)
}

# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION: pull candidate genes for pre-defined genomic regions
# Arguments:
#   region      - data frame with columns chrom, start_bp, end_bp
#   species     - Ensembl species name
#   output_file - path for the enriched output CSV
# Returns: enriched data frame (invisibly)
# ══════════════════════════════════════════════════════════════════════════════
extract_candidates_region <- function(region, species, output_file) {
  candidate_genes <- data.frame(
    chrom = character(), region_start = integer(), region_end = integer(),
    gene_id = character(), gene_start = integer(), gene_end = integer(),
    strand = integer(), biotype = character(),
    canonical_transcript = character(),
    stringsAsFactors = FALSE)

  server    <- "https://rest.ensembl.org"
  CHUNK_MAX <- 4999999L   # Ensembl REST inclusive-coord limit: end-start+1 <= 5000000

  query_chunk <- function(chrom, start, end) {
    ext <- paste0("/overlap/region/", tolower(species), "/",
                  chrom, ":", start, "-", end, "?feature=gene")
    message("  Querying: ", server, ext)
    request <- GET(paste0(server, ext), accept("application/json"))
    stop_for_status(request)
    parse_ensembl_genes(request)
  }

  for (i in seq_len(nrow(region))) {
    current_chrom <- region[i, "chrom"]
    start_bp      <- as.integer(region[i, "start_bp"])
    end_bp        <- as.integer(region[i, "end_bp"])

    # Split into ≤5 Mb chunks if the region is too large
    chunk_starts <- seq(start_bp, end_bp, by = CHUNK_MAX)
    chunk_ends   <- pmin(chunk_starts + CHUNK_MAX - 1L, end_bp)

    for (j in seq_along(chunk_starts)) {
      parsed_genes <- query_chunk(current_chrom, chunk_starts[j], chunk_ends[j])
      if (nrow(parsed_genes) > 0) {
        candidate_genes <- rbind(candidate_genes,
                                 cbind(chrom        = as.character(current_chrom),
                                       region_start = start_bp,
                                       region_end   = end_bp,
                                       parsed_genes))
      }
    }
  }
  # Remove duplicates arising from genes spanning chunk boundaries
  candidate_genes <- dplyr::distinct(candidate_genes)

  if (nrow(candidate_genes) == 0) {
    message("  No candidates found.")
    return(invisible(candidate_genes))
  }

  message("  Found ", nrow(candidate_genes), " genes — enriching via biomaRt...")
  enriched      <- enrich_with_biomart(unique(candidate_genes$gene_id), species)
  enriched_full <- dplyr::left_join(candidate_genes, enriched, by = "gene_id")
  write.table(enriched_full, file = output_file, sep = ",",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  message("  Written: ", output_file)

  invisible(enriched_full)
}

# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION: enrich a vector of Ensembl gene IDs via biomaRt
# Returns a data frame with GO terms, UniProt IDs, and Arabidopsis orthologs
# ══════════════════════════════════════════════════════════════════════════════
enrich_with_biomart <- function(gene_ids, species) {
  dataset <- if (grepl("hordeum_vulgare", tolower(species))) {
    "hvulgare_eg_gene"
  } else if (grepl("triticum_aestivum", tolower(species))) {
    "taestivum_eg_gene"
  } else {
    stop("Species not recognised — add its biomaRt dataset name manually.")
  }

  mart <- tryCatch(
    useMart("plants_mart", dataset = dataset, host = "https://plants.ensembl.org"),
    error = function(e) { warning("biomaRt connection failed: ", e$message); return(NULL) }
  )
  if (is.null(mart)) return(data.frame(gene_id = gene_ids))

  safe_getBM <- function(...) {
    tryCatch(getBM(...), error = function(e) {
      warning("getBM failed: ", e$message); return(NULL)
    })
  }

  # --- Query 1: Core annotation ---
  core <- safe_getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description",
                   "gene_biotype", "percentage_gene_gc_content"),
    filters = "ensembl_gene_id", values = gene_ids, mart = mart)

  # --- Query 2: GO terms ---
  go_raw <- safe_getBM(
    attributes = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
    filters = "ensembl_gene_id", values = gene_ids, mart = mart)

  go_collapsed <- if (!is.null(go_raw)) {
    go_raw %>%
      dplyr::filter(!is.na(go_id) & go_id != "") %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(
        go_ids     = paste(unique(go_id),          collapse = ";"),
        go_terms   = paste(unique(name_1006),      collapse = ";"),
        go_domains = paste(unique(namespace_1003), collapse = ";"),
        .groups = "drop")
  } else NULL

  # --- Query 3: Canonical transcript ---
  tx_raw <- safe_getBM(
    attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                   "transcript_is_canonical", "transcript_length", "transcript_biotype"),
    filters = "ensembl_gene_id", values = gene_ids, mart = mart)

  tx_canonical <- if (!is.null(tx_raw)) {
    tx_raw %>%
      dplyr::filter(transcript_is_canonical == 1) %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::slice_max(transcript_length, n = 1, with_ties = FALSE) %>%
      dplyr::select(ensembl_gene_id,
                    canonical_transcript_id      = ensembl_transcript_id,
                    canonical_transcript_length  = transcript_length,
                    canonical_transcript_biotype = transcript_biotype)
  } else NULL

  # --- Query 4: UniProt ---
  uniprot <- safe_getBM(
    attributes = c("ensembl_gene_id", "uniprot_gn_id"),
    filters = "ensembl_gene_id", values = gene_ids, mart = mart)

  uniprot_collapsed <- if (!is.null(uniprot)) {
    uniprot %>%
      dplyr::filter(!is.na(uniprot_gn_id) & uniprot_gn_id != "") %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(uniprot_ids = paste(unique(uniprot_gn_id), collapse = ";"),
                       .groups = "drop")
  } else NULL

  # --- Query 5: Arabidopsis orthology ---
  ara <- safe_getBM(
    attributes = c("ensembl_gene_id",
                   "athaliana_eg_homolog_ensembl_gene",
                   "athaliana_eg_homolog_associated_gene_name",
                   "athaliana_eg_homolog_orthology_type",
                   "athaliana_eg_homolog_perc_id"),
    filters = "ensembl_gene_id", values = gene_ids, mart = mart)

  ara_collapsed <- if (!is.null(ara)) {
    ara %>%
      dplyr::filter(!is.na(athaliana_eg_homolog_ensembl_gene) &
                    athaliana_eg_homolog_ensembl_gene != "") %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::slice_max(athaliana_eg_homolog_perc_id, n = 1, with_ties = FALSE) %>%
      dplyr::select(ensembl_gene_id,
                    arabidopsis_gene_id   = athaliana_eg_homolog_ensembl_gene,
                    arabidopsis_gene_name = athaliana_eg_homolog_associated_gene_name,
                    arabidopsis_orthology = athaliana_eg_homolog_orthology_type,
                    arabidopsis_perc_id   = athaliana_eg_homolog_perc_id)
  } else NULL

  # --- Join all queries ---
  result <- core %>% dplyr::rename(gene_id = ensembl_gene_id)
  if (!is.null(go_collapsed))
    result <- dplyr::left_join(result, dplyr::rename(go_collapsed,    gene_id = ensembl_gene_id), by = "gene_id")
  if (!is.null(tx_canonical))
    result <- dplyr::left_join(result, dplyr::rename(tx_canonical,    gene_id = ensembl_gene_id), by = "gene_id")
  if (!is.null(uniprot_collapsed))
    result <- dplyr::left_join(result, dplyr::rename(uniprot_collapsed, gene_id = ensembl_gene_id), by = "gene_id")
  if (!is.null(ara_collapsed))
    result <- dplyr::left_join(result, dplyr::rename(ara_collapsed,   gene_id = ensembl_gene_id), by = "gene_id")
  return(result)
}

# ══════════════════════════════════════════════════════════════════════════════
# SCRIPT BODY
# ══════════════════════════════════════════════════════════════════════════════
load("gwas_hits_master.Rdata")
ipc <- ipc %>% dplyr::arrange(desc(neglogP))

# focal_snps.Rdata should contain a character vector of SNP names (optionally
# named, e.g. c("1H" = "SNP_NAME", "2H" = "SNP_NAME2", ...)).
# Labels are taken from names() when present; otherwise derived from ipc$Chr.
load("focal_snps.Rdata")

#focal_snps<- focal_snps[5]
for (i in seq_along(focal_snps)) {
  snp_name <- focal_snps[i]
  snp_row  <- ipc %>% dplyr::filter(SNP == snp_name) %>% dplyr::select(Chr, Pos) %>% dplyr::slice(1)

  if (nrow(snp_row) == 0) {
    message("SNP '", snp_name, "' not found in ipc — skipping.")
    next
  }

  label <- if (!is.null(names(focal_snps)) && nzchar(names(focal_snps)[i])) {
    names(focal_snps)[i]
  } else {
    paste0(snp_row$Chr, "H")
  }

  out_enriched <- paste0("candidates_", label, "_enriched.csv")

  message("\n=== Locus ", label, " | ", snp_name, " | Chr", snp_row$Chr, " pos ", snp_row$Pos, " ===")
  candidates <- extract_candidates_snp(snp_row, "Hordeum_vulgare",
                                       flank_length = 2500000)

  if (nrow(candidates) == 0) {
    message("  No candidates found.")
    next
  }
  message("  Found ", nrow(candidates), " candidates — enriching via biomaRt...")
  enriched      <- enrich_with_biomart(unique(candidates$gene_id), "Hordeum_vulgare")
  enriched_full <- dplyr::left_join(candidates, enriched, by = "gene_id")
  write.table(enriched_full, file = out_enriched, sep = ",",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  message("  Written: ", out_enriched)
}














# ── Region query example: 6H 545–550 Mb ──────────────────────────────────────
test_region <- data.frame(
  chrom    = "6H",
  start_bp = 541.4e6,
  end_bp   = 559.2e6
)

region_results <- extract_candidates_region(
  test_region, "Hordeum_vulgare",
  output_file = "candidates_region_6H_enriched_wide.csv"
)

# ── Region query example: 5H 22.6–27.6 Mb ──────────────────────────────────────
test_region <- data.frame(
  chrom    = "5H",
  start_bp = 22.6e6,
  end_bp   = 27.6e6
)

region_results <- extract_candidates_region(
  test_region, "Hordeum_vulgare",
  output_file = "candidates_region_5H_enriched1.csv"
)

#  1 JHI-Hv50k-2016-286560 5

# ── Region query example: 5H 14.9–19.6 Mb ──────────────────────────────────────
test_region <- data.frame(
  chrom    = "5H",
  start_bp = 14.9e6,
  end_bp   = 19.6e6
)

region_results <- extract_candidates_region(
  test_region, "Hordeum_vulgare",
  output_file = "candidates_region_5H_enriched2.csv"
)

#  1 JHI-Hv50k-2016-283964 5
