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
    cue = tar_cue(mode = "thorough"),
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
  
  parse_BU <- function(counts_path, biospecimen, cell_index_df) {
    # read counts
    raw_counts <- read.csv(counts_path, row.names = 1, check.names = FALSE) |> as.matrix()
    
    # fix names
    colnames(raw_counts) <- gsub("\\.([0-9]+)$", "-\\1", colnames(raw_counts)) # restore 10x barcode format
    rownames(raw_counts) <- gsub("\\.\\d+$", "", rownames(raw_counts))         # remove gene version suffix
    
    # create SCE
    sce = SingleCellExperiment(assays = list(counts = raw_counts))
    sce$sample_id <- biospecimen
    
    cd <- colData(sce) |> 
      as.data.frame() |> 
      tibble::rownames_to_column(var = "old_cell_id")
    
    # inner_join to filter
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
    
    # 2. Read processed files to get filtered cell IDs
    tar_target(
      bu_processed_df,
      {
        set.seed(19980)
        file_metadata |> 
          filter(Atlas.Name=="HTAN BU") |>
          filter(Biospecimen |> str_detect(",")) |> 
          pull(full_path) |> 
          map_dfr(~ read.csv(.x) |> dplyr::slice(-1)) |> 
          
        # Create cell Id index
          group_by(SampleID) |>
          mutate(cell_index = row_number()) |>
          ungroup()
        
        
      },
      packages = c("dplyr", "purrr", "readr")
    ),
    
    # 3. Get filtered cell IDs
    tar_target(
      cell_index_df,
      {
        bu_processed_df |> 
          # This is ridicular to understand. For some reason, SampleID missed a digital from Biospecimen. 
          # i.e, SampleID HTA3_8001_001, Biospecimen HTA3_8001_1001
          mutate(
            SampleID = sub("_(\\d+)$", "_1\\1", SampleID)
          ) |>
          select(NAME, cell_index, SampleID)
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
          cell_index_df = cell_index_df
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
          "/vast/scratch/users/shen.m/htan/hta_2025/0.1.0/counts/"
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

# One anndata
library(zellkonverter)
anndata = readH5AD("/vast/scratch/users/shen.m/htan/hta_2025/0.1.0/counts/HTA3_8001_1001.h5ad",reader = "R", use_hdf5 = T)
anndata
