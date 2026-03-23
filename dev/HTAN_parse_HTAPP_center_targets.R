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
  
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue",
    workspace_on_error = TRUE,
    cue = tar_cue(mode = "never"),
    controller = crew_controller_slurm(
      name = "elastic",
      workers = 300,
      tasks_max = 20,
      seconds_idle = 30,
      crashes_error = 10,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = c(30, 40, 60, 100, 150), 
        cpus_per_task = c(2),
        time_minutes = c(60*24),
        verbose = T
      )
    )
  )
  
  # ------------------------------------------------------------
  # Helper functions
  # ------------------------------------------------------------
  
  parse_HTAPP <- function(sample_id, mtx_path, genes_path, barcodes_path, processed_path) {
    counts <- ReadMtx(mtx = mtx_path, 
                      cells = barcodes_path, 
                      features = genes_path, 
                      # Use Ensemble 
                      feature.column = 1)
    
    genes <- read.delim(genes_path, header = FALSE)$V1 # Extract Ensemble ID 
    barcodes <- read.delim(barcodes_path, header = FALSE)$V1 
    filtered_df <- read.table(processed_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)[-1,]
    filtered_cell_id <- filtered_df |> pull(NAME)
    
    
    sce <- SingleCellExperiment(
      assays = list(counts = counts),
      rowData = data.frame(gene_id = genes) |> tibble::column_to_rownames("gene_id"),
      colData = data.frame(barcode = barcodes,
                           sample_id = sample_id)
    )
    
    sce_subset <- sce[, colnames(sce) %in% filtered_cell_id]
    
    colnames(sce_subset) <- paste0(sce_subset$barcode, "___", sample_id)
    
    sce_subset
  }
  
  save_h5ad <- function(sample_id, sce, save_directory) {
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
    
    # 3. Group by sample_id for parallel processing
    tar_target(
      htapp_metadata_grouped,
      htapp_metadata |>
        group_by(sample_id) |>
        tar_group(),
      iteration = "group"
    ),
    
    # 4. Process each sample in parallel
    tar_target(
      htapp_sce_per_sample,
      {
        row <- htapp_metadata_grouped
        # After pivot_wider, each row should have all file paths
        parse_HTAPP(
          sample_id = row$sample_id[1],
          mtx_path = row$mtx_path[1],
          genes_path = row$genes_path[1],
          barcodes_path = row$barcodes_path[1],
          processed_path = row$processed_file[1]
        )
      },
      pattern = map(htapp_metadata_grouped),
      iteration = "list",
      packages = c("SingleCellExperiment", "Seurat", "dplyr", "tibble")
    ),
    
    # 5. Create tibble with sample_id and sce
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
    
    # 6. Save to disk
    tar_target(
      save_htapp_h5ad,
      {
        if (is.null(htapp_sce_tbl$sce[[1]])) return(NULL)
        save_h5ad(
          htapp_sce_tbl$sample_id, 
          htapp_sce_tbl$sce[[1]], 
          "~/scratch/cache_temp/hta/09-11-2025/counts/"
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

tar_workspace(save_htapp_h5ad_b9d77b97beed2a7a, store = store, script = paste0(store, "_target_script.R") )
tar_meta(starts_with("sce_list_by_sample_list"), store = store) |> filter(is.na(error)) |> pull(name)
sce_list_by_sample_list = tar_read(sce_list_by_sample_list, store = store) |> bind_rows()


