library(targets)
library(dplyr)
store = "~/scratch/htan/BU_target_store"

tar_script({
  
  library(targets)
  library(tarchetypes)
  library(tidyverse)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(tidyr)
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
  
  parse_BU <- function(counts_path, biospecimen, filtered_cell_ids) {
    # read counts
    raw_counts <- read.csv(counts_path, row.names = 1, check.names = FALSE) |> as.matrix()
    
    # fix names
    colnames(raw_counts) <- gsub("\\.([0-9]+)$", "-\\1", colnames(raw_counts)) # restore 10x barcode format
    rownames(raw_counts) <- gsub("\\.\\d+$", "", rownames(raw_counts))         # remove gene version suffix
    
    # create SCE
    sce = SingleCellExperiment(assays = list(counts = raw_counts))
    sce$sample_id <- biospecimen
    
    # Filter cells
    sce_subset <- sce[, colnames(sce) %in% filtered_cell_ids]
    
    # Update colnames with sample_id
    colnames(sce_subset) <- paste0(colnames(sce_subset), "___", biospecimen)
    
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
    
    # 2. Read processed files to get filtered cell IDs
    tar_target(
      bu_processed_df,
      {
        file_metadata |> 
          filter(Atlas.Name=="HTAN BU") |>
          filter(Biospecimen |> str_detect(",")) |> 
          pull(full_path) |> 
          map_dfr(~ read.csv(.x) |> dplyr::slice(-1))
      },
      packages = c("dplyr", "purrr", "readr")
    ),
    
    # 3. Get filtered cell IDs
    tar_target(
      bu_filtered_cell_ids,
      {
        bu_processed_df |> pull(NAME)
      }
    ),
    
    # 4. Process BU metadata for raw counts files
    tar_target(
      bu_metadata,
      {
        file_metadata |> 
          filter(Atlas.Name=="HTAN BU") |> 
          filter(full_path |> str_detect("raw")) |>
          dplyr::select(full_path, Biospecimen) |>
          mutate(sample_id = Biospecimen)
      },
      packages = c("dplyr", "stringr")
    ),
    
    # 5. Group by sample_id for parallel processing
    tar_target(
      bu_metadata_grouped,
      bu_metadata |>
        group_by(sample_id) |>
        tar_group(),
      iteration = "group"
    ),
    
    # 6. Process each sample in parallel
    tar_target(
      bu_sce_per_sample,
      {
        row <- bu_metadata_grouped
        # Each group may have multiple files, process each file
        # If multiple files per sample, we might need to merge them
        # For now, assuming one file per sample_id
        parse_BU(
          counts_path = row$full_path[1],
          biospecimen = row$sample_id[1],
          filtered_cell_ids = bu_filtered_cell_ids
        )
      },
      pattern = map(bu_metadata_grouped),
      iteration = "list",
      packages = c("SingleCellExperiment", "dplyr")
    ),
    
    # 7. Create tibble with sample_id and sce
    tar_target(
      bu_sce_tbl,
      {
        if (is.null(bu_sce_per_sample)) return(NULL)
        
        tibble(
          sample_id = unique(colData(bu_sce_per_sample)$sample_id),
          sce = list(bu_sce_per_sample)
        )
      },
      pattern = map(bu_sce_per_sample),
      iteration = "list"
    ),
    
    # 8. Save to disk
    tar_target(
      save_bu_h5ad,
      {
        if (is.null(bu_sce_tbl$sce[[1]])) return(NULL)
        save_h5ad(
          bu_sce_tbl$sample_id, 
          bu_sce_tbl$sce[[1]], 
          "~/scratch/cache_temp/hta/09-11-2025/counts/"
        )
      },
      pattern = map(bu_sce_tbl),
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



