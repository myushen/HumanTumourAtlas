# =============================================================================
# Parse HTAN Data Files to SingleCellExperiment Objects
# =============================================================================
#
# This script parses single-cell RNA sequencing data for lung tissues from the Human Tumor Atlas
# Network (HTAN) into SCE objects. It handles two data
# formats:
#
# 1. CSV format (uncompressed):
#    - Files: *.mtx, *_counts_genes.csv, *_counts_barcodes.csv
#    - Gene identifiers: Gene symbols (converted to Ensembl IDs)
#    - Processing: Maps gene symbols to Ensembl IDs using EnsDb.Hsapiens.v86
#
# 2. TSV format (compressed):
#    - Files: *_matrix.mtx.gz, *_features.tsv.gz, *_barcodes.tsv.gz
#    - Gene identifiers: Ensembl IDs (already present)
#    - Processing: Direct use of Ensembl IDs, no conversion needed
#
# Key features:
# - Handles duplicate genes and barcodes by keeping first occurrence
# - Creates unique cell IDs: barcode___Biospecimen
# - Supports multi-channel data (combines channels per biospecimen)
# - Filters out Level 4 (cell type, H5AD, cluster etc) processed data from HTAN metadata
#
# Output: SingleCellExperiment objects with:
#   - assays$counts: Gene expression count matrix (Ensembl IDs as rownames)
#   - colData: Cell metadata including barcode, Biospecimen, and cell_id
#   - rowData: Gene metadata with Ensembl gene IDs
# =============================================================================
library(ensembldb) 
library(dplyr)
library(purrr)
library(stringr)
library(fs)
library(readr)
library(Matrix)
library(R.utils)
library(tidyr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(zellkonverter)
library(Seurat)
library(tibble)

parse_csv <- function(sample, mtx_path, genes_path, barcodes_path) {
  # --- Load input files ---
  counts <- readMM(mtx_path) |> t()
  genes <- read.csv(genes_path, header = FALSE)$V2
  barcodes <- read.csv(barcodes_path, header = FALSE)$V2
  
  stopifnot(length(genes) == nrow(counts))
  stopifnot(length(barcodes) == ncol(counts))
  
  # --- Deduplicate genes ---
  dup_genes <- duplicated(genes)
  if (any(dup_genes)) {
    keep_idx <- !duplicated(genes)
    counts <- counts[keep_idx, , drop = FALSE]
    genes <- genes[keep_idx]
  } 
  
  # --- Deduplicate barcodes ---
  dup_barcodes <- duplicated(barcodes)
  if (any(dup_barcodes)) {
    keep_idx <- !duplicated(barcodes)
    counts <- counts[, keep_idx, drop = FALSE]
    barcodes <- barcodes[keep_idx]
  } 
  
  # --- Create SCE object ---
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    rowData = data.frame(gene_id = genes) |> tibble::column_to_rownames("gene_id"),
    colData = data.frame(barcode = barcodes,
                         #counts_path = basename(mtx_path),
                         sample_id = sample)
    
  )
  
  colnames(sce) <- paste0(sce$barcode, "___", sample)
  
  # --- Convert SYMBOL to Ensembl ID
  gene_id <- ensembldb::mapIds(
    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    keys       = rownames(sce), 
    column     = "GENEID",      # output: Ensembl gene id
    keytype    = "GENENAME",    # input: gene symbol
    multiVals  = "first"        # choose first mapping when duplicates
  )
  
  # drop unmapped
  keep <- !is.na(gene_id)
  sce <- sce[keep, ]
  gene_id <- gene_id[keep]
  rowData(sce)$gene_id <- gene_id
  rownames(sce) <- gene_id
  
  sce
}


parse_tsv <- function(sample, mtx_path, genes_path, barcodes_path) {
  counts <- ReadMtx(mtx = mtx_path, 
                    cells = barcodes_path, 
                    features = genes_path, 
                    # Use Ensemble 
                    feature.column = 1)
  
  genes <- read.delim(genes_path, header = FALSE)$V1 # Extract Ensemble ID 
  barcodes <- read.delim(barcodes_path, header = FALSE)$V1 
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    rowData = data.frame(gene_id = genes) |> tibble::column_to_rownames("gene_id"),
    colData = data.frame(barcode = barcodes,
                         #counts_path = basename(mtx_path),
                         sample_id = sample)
    
  )
  
  colnames(sce) <- paste0(sce$barcode, "___", sample)
  
  sce
    
}

metadata_path = "/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv"
file_path = "/vast/scratch/users/shen.m/synapse_data/lung/counts/"

csv_files <- read.csv(metadata_path,sep = "\t", na.strings = c("NA",""), header = TRUE) |>
  filter(Level != "Level 4", # exclude processed
         str_detect(Filename, "mtx$|genes.csv$|barcodes.csv$")) |>   
  mutate(
    Filename_basename = basename(Filename),
    full_path = file.path(file_path, Filename_basename)
  ) |> mutate(
    file_type = case_when(
      str_detect(Filename_basename, "molecule_counts\\.mtx$") ~ "mtx_path",
      str_detect(Filename_basename, "_counts_genes\\.csv$") ~ "genes_path",
      str_detect(Filename_basename, "_counts_barcodes\\.csv$") ~ "barcodes_path",
      TRUE ~ NA_character_
    )
  ) |>
  select(Biospecimen, file_type, full_path, everything()) |>
  pivot_wider(
    id_cols = Biospecimen,
    names_from = file_type,
    values_from = full_path
  ) |> 
  mutate(sample_id = Biospecimen) |> 
  
  # ONE SAMPLE FOR TESTING PURPOSE
  head(1)

tsv_files <- read.csv(metadata_path,sep = "\t", na.strings = c("NA",""), header = TRUE) |>
  filter(Level != "Level 4",  # exclude processed
         str_detect(Filename, "mtx.gz$|barcodes.tsv.gz$|features.tsv.gz$")) |>  
  mutate(
    Filename_basename = basename(Filename),
    full_path = file.path(file_path, Filename_basename)
  ) |> mutate(
    file_type = case_when(
      str_detect(Filename_basename, "_matrix.mtx.gz$") ~ "mtx_path",
      str_detect(Filename_basename, "_features.tsv.gz$") ~ "genes_path",
      str_detect(Filename_basename, "_barcodes.tsv.gz$") ~ "barcodes_path",
      TRUE ~ NA_character_
    ),
    #channel_key = str_extract(Filename_basename, "^.*_channel[0-9]+"),
    channel_number = str_extract(Filename_basename, "channel[0-9]+")
  ) |>
  select(Biospecimen,channel_number, file_type, full_path, everything()) |>
  pivot_wider(
    id_cols = c(Biospecimen,channel_number),
    names_from = file_type,
    values_from = full_path
  ) |> 
  mutate(sample_id = paste(Biospecimen, channel_number, sep = "___")) |> 
  
  # ONE SAMPLE FOR TESTING PURPOSE
  head(1)
  

# =============================================================================
# Create SCE object
# =============================================================================
save_directory = "/vast/scratch/users/shen.m/htan/"

# From a uncompressed csv input
parse_csv(csv_files$Biospecimen, csv_files$mtx_path, csv_files$genes_path, csv_files$barcodes_path) |> 
  writeH5AD(file.path(save_directory, paste0(csv_files$sample_id, ".h5ad")), compression = "gzip")

# From a compressed tsv input
parse_tsv(tsv_files$sample_id, tsv_files$mtx_path, tsv_files$genes_path, tsv_files$barcodes_path) |> 
  writeH5AD(file.path(save_directory, paste0(tsv_files$sample_id, ".h5ad")), compression = "gzip")

# sce2 = tsv_files |> mutate(sce = pmap(
#   list(
#     tsv_files$sample_id, tsv_files$mtx_path, tsv_files$genes_path, tsv_files$barcodes_path
#   ), parse_tsv)
# ) |> 
#   group_by(Biospecimen) |>
#   summarise(sce_list = list(sce), .groups = "drop") |> 
#   mutate(sce = map(sce_list, function(x) {
#     x <- unlist(x, recursive = T)
#     if (length(x) == 1) return(x[[1]])
#     Reduce(cbind, x)
#   }))

# sce2 |> 
#   pwalk(~ writeH5AD(
#     ..3, 
#     file.path(save_directory, paste0(..1, ".h5ad")),
#     compression = "gzip"
#     )
#   )

# # =============================================================================
# # Check saved SCE object
# # =============================================================================
# x = readH5AD(file.path(save_directory, "HTA1_203_332101___channel3.h5ad"), reader = "R")
# y = readH5AD(file.path(save_directory, "HTA8_4003_001101.h5ad"), reader = "R")


           