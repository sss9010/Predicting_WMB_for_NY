# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This directory holds raw and processed genotyping output for the **OREI (Organic Rye Evaluation Initiative) 2024** barley genotyping batches. Data were generated on the **Barley 3K (B3K) Illumina SNP array**, processed through **Illumina GenomeStudio / BrainScape**, and variants are mapped to the **Morex V3** barley reference genome.

This data feeds into the downstream genomic prediction analyses in the parent `NY_WMBCL` project — see the parent `CLAUDE.md` for that pipeline context.

## Batch Structure

Three independent genotyping batches, each with the same file layout:

| Batch | Samples |
|-------|---------|
| `OREI2024-B3K_01-10_Tri/` | Plates 01–10, Tri experiment |
| `OREI2024-B3K_11-14_DualO/` | Plates 11–14, DualO experiment |
| `OREI2024-B3K_01-14-COMBO/` | Combined plates 01–14 |

## File Formats per Batch

| File | Format | Contents |
|------|--------|----------|
| `*-AB.txt` | Tab-delimited | Allele calls in AB format (columns: ID, CHROM, POS, then ~960 sample genotypes as AA/AB/BB/NC) |
| `*-FDT.txt` | Tab-delimited | Full Detailed Table — per-sample columns: Genotype, GC Score, Theta, R |
| `*-mFDT.txt` | Tab-delimited | Minimal FDT — per-sample columns: GT, Theta, R only |
| `*-R.txt` | Tab-delimited | Raw R intensity values (signal strength) |
| `*-Samples.txt` | Tab-delimited | Sample metadata: Index, Well, Sample ID, Plate, Call Rate, GC Score, Sentrix ID/Position |
| `*-VCF-Samples.txt` | Plain text | Sample IDs passing QC |
| `*-poorQC_VCF-Samples.txt` | Plain text | Sample IDs failing QC |
| `*.MorexV3.vcf` | VCF 4.x | Variant calls mapped to Morex V3 |
| `*.MorexV3.IDnum.vcf` | VCF 4.x | Same as above with numeric variant IDs |
| `*-Project/*.bsc` | Binary | GenomeStudio project file |
| `*-Project/Data/*.bin` | Binary | GenomeStudio internal data (ld, ad, heredity, clusters, etc.) — not directly parseable |

## Key Notes

- **Chromosomes**: Barley diploid (2n=14), chromosomes named `1H`–`7H` in CHROM field.
- **Sample IDs**: Follow `BAROREI[number]` convention.
- **QC threshold**: Samples are split into passing (`*-VCF-Samples.txt`) and failing (`*-poorQC_VCF-Samples.txt`) based on Call Rate in `*-Samples.txt`.
- **COMBO batch**: Contains the merged dataset across all plates; use this for full-cohort analyses rather than combining Tri + DualO manually.
- **`.bin` files**: Internal to GenomeStudio/BrainScape and not meant for direct parsing — use the exported `.txt` or `.vcf` files for downstream R analysis.

## Downstream Integration

For genomic prediction, the `.vcf` files are the primary input to `ASRgenomics` / `rrBLUP` pipelines in the parent project's `data/geno/` directory. The standard path is:

1. Filter samples using `*-VCF-Samples.txt`
2. Use `*.MorexV3.vcf` (or `.IDnum.vcf` for numeric marker IDs) as genotype input
3. Convert to HapMap or numeric matrix format per the parent pipeline conventions
