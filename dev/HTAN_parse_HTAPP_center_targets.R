library(targets)
library(dplyr)
store = "~/scratch/htan/HTAPP_target_store"

tar_script({
  
  library(targets)
  library(tarchetypes)
  library(tidyverse)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(tidyr)
  library(Seurat)
  library(crew)
  library(crew.cluster)
  
  # ------------------------------------------------------------
  # Target options
  # ------------------------------------------------------------
  
  # Helper (optional) to avoid repetition
  new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 2, backup = NULL) {
    crew_controller_slurm(
      name = name,
      workers = workers,
      crashes_max = crashes_max,
      seconds_idle = 30,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = mem_gb,
        cpus_per_task = cpus_per_task,
        time_minutes = time_min
      ),
      backup = backup
    )
  }
  
  # Small → large, with fallbacks to the next size up
  elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2)
  elastic_120  <- new_elastic("elastic_120",  120,  60 * 4,  workers = 16, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
  elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
  elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
  elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
  elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
  
  elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)
  
  # Group for targets (small → large)
  controllers <- crew_controller_group(
    elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_5_minimal
  )
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    cue = tar_cue(mode = "thorough"),
    format = "qs",
    #debug = "dataset_id_sct_ea377f6e2d0ae2b7",
    workspace_on_error = TRUE,
    controller = controllers, 
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    ) 
  )
  
  # ------------------------------------------------------------
  # Helper functions
  # ------------------------------------------------------------
  
  parse_HTAPP <- function(sample_id, mtx_path, genes_path, barcodes_path, cell_index_df) {
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
                           sample_id = sample_id)
    )
    
    cd <- colData(sce) |> 
      as.data.frame() |> 
      tibble::rownames_to_column(var = "old_cell_id")
    
    cd <- cd |> inner_join(cell_index_df, by = c("old_cell_id" = "NAME")) |> 
      dplyr::rename(cell_id = cell_index)
    
    sce_subset = sce[, colnames(sce) %in% cell_index_df$NAME]
    stopifnot(nrow(cd) == ncol(sce_subset))
    
    colnames(sce_subset) <- as.numeric(cd$cell_id)
    
    sce_subset <- if (all(grepl("^ENSG", rownames(sce_subset)))) {
      sce_subset
    } else {
      sce_subset |> convert_gene_to_ensemble()
    }
    sce_subset
  }
  
  convert_gene_to_ensemble <- function(sce) {
    # --- Convert SYMBOL to Ensembl ID
    library(AnnotationFilter)
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
  
  save_h5ad <- function(sample_id, sce, save_directory) {
    if (!dir.exists(save_directory)) dir.create(save_directory, recursive = T)
    zellkonverter::writeH5AD(sce, file = paste0(save_directory, sample_id, ".h5ad"),
                             compression = "gzip")
    print(paste("saved successfully:", sample_id))
  }
  
  # ------------------------------------------------------------
  # Pipeline
  # ------------------------------------------------------------
  list(
    
    # 1. Read file metadata
    tar_target(
      file_metadata,
      {
        file_path = "/vast/scratch/users/shen.m/synapse_data/lung/counts/"
        
        file_metadata <- read.csv("/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv",
                                  sep = "\t", na.strings = c("NA",""), header = TRUE) |> as_tibble() |> 
          # THIS IS ASSIGNED WRONG TO BIOSPECIMEN. HTAN PHASE1 IS SOOOO COMPLEX!
          dplyr::filter(Filename != "single_cell_RNAseq_level_4_lung/lung_HTA1_203_332102_ch1_L4.tsv") |>
          mutate(
            Filename_basename = basename(Filename),
            full_path = file.path(file_path, Filename_basename))
        
        file_metadata
      }
    ),
    
    # 2. Process HTAPP metadata
    tar_target(
      htapp_metadata,
      {
        file_metadata |> 
          filter(Atlas.Name=="HTAN HTAPP") |> 
          arrange(desc(Biospecimen)) |> 
          mutate(
            file_type = case_when(
              str_detect(Filename_basename, "_matrix.mtx.gz$") ~ "mtx_path",
              str_detect(Filename_basename, "_features.tsv.gz$") ~ "genes_path",
              str_detect(Filename_basename, "_barcodes.tsv.gz$") ~ "barcodes_path",
              str_detect(Filename_basename, "_L4.tsv$") ~ "processed_file"
            ),
            channel_number = Filename_basename |>
              str_replace_all("ch(?=[0-9]+)", "channel") |>
              str_extract("channel[0-9]+")
          ) |> 
          group_split(Biospecimen) |>
          map_dfr(~ {
            subgroup <- .x
            if (nrow(subgroup) > 4) {
              subgroup <- subgroup |>
                mutate(sample_id = paste(Biospecimen, channel_number, sep = "___"))
            } else {
              subgroup <- subgroup |>
                mutate(sample_id = Biospecimen)
            }
            subgroup
          }) |>
          dplyr::select(full_path, Biospecimen, sample_id, file_type) |>
          pivot_wider(
            id_cols = c(Biospecimen, sample_id),
            names_from = file_type,
            values_from = full_path
          )
      },
      packages = c("dplyr", "tidyr", "stringr", "purrr")
    ),
    
    # 3. Build global NAME -> cell_index map keyed by sample_id
    tar_target(
      htapp_cell_index_df,
      {
        htapp_metadata |>
          dplyr::select(sample_id, processed_file) |>
          dplyr::distinct() |>
          purrr::pmap_dfr(function(sample_id, processed_file) {
            read.table(
              processed_file,
              sep = "\t",
              header = TRUE,
              stringsAsFactors = FALSE
            )[-1, ] |>
              dplyr::mutate(
                sample_id = sample_id,
                cell_index = dplyr::row_number()
              ) |>
              dplyr::select(sample_id, NAME, cell_index)
          })
      },
      packages = c("dplyr", "purrr")
    ),
    
    # 4. Group by sample_id for parallel processing
    tar_target(
      htapp_metadata_grouped,
      htapp_metadata |>
        group_by(sample_id) |>
        tar_group(),
      iteration = "group"
    ),
    
    # 5. Process each sample in parallel
    tar_target(
      htapp_sce_per_sample,
      {
        row <- htapp_metadata_grouped
        sample_cell_index_df <- htapp_cell_index_df |>
          dplyr::filter(sample_id == row$sample_id[1]) |>
          dplyr::select(NAME, cell_index)
        parse_HTAPP(
          sample_id = row$sample_id[1],
          mtx_path = row$mtx_path[1],
          genes_path = row$genes_path[1],
          barcodes_path = row$barcodes_path[1],
          cell_index_df = sample_cell_index_df
        )
      },
      pattern = map(htapp_metadata_grouped),
      iteration = "list",
      packages = c("SingleCellExperiment", "Seurat", "dplyr", "tibble")
    ),
    
    # 6. Create tibble with sample_id and sce
    tar_target(
      htapp_sce_tbl,
      {
        if (is.null(htapp_sce_per_sample)) return(NULL)
        
        tibble(
          sample_id = unique(colData(htapp_sce_per_sample)$sample_id),
          sce = list(htapp_sce_per_sample)
        )
      },
      pattern = map(htapp_sce_per_sample),
      iteration = "list"
    ),
    
    # 7. Save to disk
    tar_target(
      save_htapp_h5ad,
      {
        if (is.null(htapp_sce_tbl$sce[[1]])) return(NULL)
        save_h5ad(
          htapp_sce_tbl$sample_id, 
          htapp_sce_tbl$sce[[1]], 
          "/vast/scratch/users/shen.m/htan/hta_2025/0.1.0/counts/"
        )
      },
      pattern = map(htapp_sce_tbl),
      packages = c("zellkonverter")
    )
  )
  
}, script = paste0(store, "_target_script.R"), ask = FALSE)


job::job({
  
  tar_make(
    script = paste0(store, "_target_script.R"), 
    store = store, 
    reporter = "summary"
  )
  
})

# One anndata
library(zellkonverter)
anndata = readH5AD("/vast/scratch/users/shen.m/htan/hta_2025/0.1.0/counts/HTA1_203_332101.h5ad",reader = "R", use_hdf5 = T)
anndata

