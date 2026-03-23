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
library(targets)
library(dplyr)
store_file_hta_lung = "/vast/scratch/users/shen.m/htan/hta_lung_scrnaseq_target_store"

tar_script({
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
  library(targets)
  library(crew)
  library(crew.cluster)
  
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue",
    workspace_on_error = TRUE,
    cue = tar_cue(mode = "thorough"),
    controller = crew_controller_slurm(
      name = "elastic",
      workers = 300,
      tasks_max = 20,
      seconds_idle = 30,
      crashes_error = 10,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = c(20, 40, 60, 100, 150), 
        cpus_per_task = c(2),
        time_minutes = c(60*24),
        verbose = T
      )
    )
  )
  
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
  
  read_metadata_file <- function(metadata_path, file_path, format ){
    if (format == "csv") {
      read.csv(metadata_path,sep = "\t", na.strings = c("NA",""), header = TRUE) |>
        dplyr::filter(Level != "Level 4", # exclude processed
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
        dplyr::select(Biospecimen, file_type, full_path, everything()) |>
        pivot_wider(
          id_cols = Biospecimen,
          names_from = file_type,
          values_from = full_path
        ) |> 
        mutate(sample_id = Biospecimen)
    } else if (format == "tsv") {
      read.csv(metadata_path,sep = "\t", na.strings = c("NA",""), header = TRUE) |>
        dplyr::filter(Level != "Level 4",  # exclude processed
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
        dplyr::select(Biospecimen,channel_number, file_type, full_path, everything()) |>
        pivot_wider(
          id_cols = c(Biospecimen,channel_number),
          names_from = file_type,
          values_from = full_path
        ) |> 
        mutate(sample_id = paste(Biospecimen, channel_number, sep = "___"))
    }
  }
  
  save_h5ad_to_directory <- function(files, save_directory, format) {
    # When used with pattern = map(), files is a single row (tibble/data.frame with 1 row)
    # Extract values as scalars and return the file path for format = "file"
    
    output_file <- file.path(save_directory, paste0(files$sample_id[[1]], ".h5ad"))
    
    if (format == "csv") {
      parse_csv(
        sample = files$sample_id[[1]], 
        mtx_path = files$mtx_path[[1]], 
        genes_path = files$genes_path[[1]], 
        barcodes_path = files$barcodes_path[[1]]
      ) |> 
        writeH5AD(
          output_file, 
          compression = "gzip"
        )
      
      print("H5AD saved successfully .. ")
    } else if (format == "tsv") {
      parse_tsv(
        sample = files$sample_id[[1]], 
        mtx_path = files$mtx_path[[1]], 
        genes_path = files$genes_path[[1]], 
        barcodes_path = files$barcodes_path[[1]]
      ) |> 
        writeH5AD(
          output_file, 
          compression = "gzip"
        )
      
      print("H5AD saved successfully .. ")
    }
    
    # Return the file path - required when using format = "file"
    return(output_file)
  }
  
  list(
    tar_target(metadata_path, "/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv"),
    tar_target(file_path, "/vast/scratch/users/shen.m/synapse_data/lung/counts/"),
    tar_target(
      csv_files,
      read_metadata_file(metadata_path, file_path, format = "csv") 
        # # FOR TESTING PURPOSE
        # head(2),
    ),
    tar_target(
      tsv_files,
      read_metadata_file(metadata_path, file_path, format = "tsv") 
        # # FOR TESTING PURPOSE
        # head(2)
    ),
    tar_target(
      saved_h5ad_from_csv,
      save_h5ad_to_directory(csv_files, "/vast/scratch/users/shen.m/hta/09-11-2025/counts", "csv" ),
      pattern = map(csv_files),
      format = "file"
    ),
    
    tar_target(
      saved_h5ad_from_tsv,
      save_h5ad_to_directory(tsv_files, "/vast/scratch/users/shen.m/hta/09-11-2025/counts", "tsv" ),
      pattern = map(tsv_files),
      format = "file"
    )
  )
  
}, script = paste0(store_file_hta_lung, "_target_script.R"), ask = FALSE)


job::job({
  
  tar_make(
    script = paste0(store_file_hta_lung, "_target_script.R"), 
    store = store_file_hta_lung, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

# tar_workspace(saved_h5ad_from_csv_7cea16f9edd8e615,store = store_file_hta_lung, script = paste0(store_file_hta_lung, "_target_script.R") )
# tar_meta(store = store_file_hta_lung) |> dplyr::filter(!is.na(error)) |> distinct(name, error)
# debugonce(save_h5ad_to_directory)
# save_h5ad_to_directory(tsv_files, "/vast/scratch/users/shen.m/hta/09-11-2025/counts", "tsv" )
# tsv_files <- tar_read(tsv_files,store = store_file_hta_lung )
