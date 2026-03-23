library(targets)
library(dplyr)
store = "~/scratch/htan/MSK_target_store"

tar_script({
  
  library(targets)
  library(tarchetypes)
  library(tidyverse)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(tidyr)
  library(tidySingleCellExperiment)
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
  # process one file
  # ------------------------------------------------------------
  process_one_file <- function(fp, sample_id) {
    
    map_id <- read.csv("/vast/scratch/users/shen.m/synapse_data/lung/counts/adata_sample_id_htan_id_map.csv")
    
    
    data <- zellkonverter::readH5AD(fp, reader = "R", use_hdf5 = TRUE)
    
    if (!"patient" %in% colnames(as.data.frame(colData(data))))
      return(NULL)
    
    # Lowest level, because other columns will be append in cellNexus API
    wanted <-  "patient"
    # |> grep("patient|cell_type", colnames(colData(sce)), value = TRUE)
    
    colData(data) <- colData(data)[ , wanted, drop = FALSE]
    
    # extract sample + barcode
    data <- data %>%
      tidyr::separate(".cell",
                      c("sample", "barcode"),
                      sep = "_(?=[^_]+$)",
                      remove = FALSE) %>%
      left_join(map_id, by = "sample") %>%
      dplyr::rename(sample_id = sample_HTAN_ID)
    
    # ---- safe slicing so NA does not break ----
    keep <- which(colData(data)$sample_id %in% sample_id)
    if (length(keep) == 0) return(NULL)
    data <- data[, keep]
    
    # keep first assay only
    assay_name <- names(assays(data))[1]
    assays(data) <- assays(data)[assay_name]
    assayNames(data) <- "counts"
    
    # Remove extra layers
    reducedDims(data) <- NULL
    metadata(data) <- list()
    
    data
  }
  

  save_h5ad <- function(sample_id, sce, save_directory) {
    .x = sce[[1]]
    zellkonverter::writeH5AD(.x, file = paste0(save_directory, sample_id, ".h5ad"),
                             compression = "gzip")
    print("saved successfully..")
  }
  
  # ------------------------------------------------------------
  # Pipeline
  # ------------------------------------------------------------
  list(
    
    # 1. Read metadata
    tar_target(
      h5ad_metadata,
      {
        h5ad_list <- list.files("/vast/scratch/users/shen.m/synapse_data/lung/counts/", pattern = ".h5ad", full.name = TRUE) |> 
          map(~ {
            data = zellkonverter::readH5AD(.x, reader= "R", use_hdf5 = T)
            data$full_path = .x
            data
          }, .progress = T)
        
        map_id <- read.csv("/vast/scratch/users/shen.m/synapse_data/lung/counts/adata_sample_id_htan_id_map.csv")
        
        h5ad_metadata <- h5ad_list |> map( ~ colData(.x) |> as.data.frame() |> rownames_to_column(".cell"), .progress = T) |> 
          bind_rows() |> as_tibble() |> tidyr::separate(".cell", c("sample", "barcode"), sep = "_(?=[^_]+$)",  remove=FALSE) |>
          left_join(map_id, by = "sample", copy=TRUE) |> 
          dplyr::rename(sample_id = sample_HTAN_ID) |> 
          filter(!is.na(sample_id)) |> 
          separate_rows(sample_id, sep = ";") |>
          distinct(sample_id, full_path) |> 
          add_count(sample_id)
        
        h5ad_metadata
      }
    ),
    
    tar_target(
      metadata_group_by_sample_id,
      h5ad_metadata     |>
        group_by(sample_id, full_path) |>
        tar_group(),
      iteration = "group"
    ),
    
    # 4. Process all files for this sample
    tar_target(
      sce_list_per_sample,
      process_one_file(metadata_group_by_sample_id$full_path, metadata_group_by_sample_id$sample_id),
      pattern = map(metadata_group_by_sample_id),
      iteration = "list"
    ),
    
    tar_target(
      sce_list_by_sample_list,
      {
        
        if (sce_list_per_sample |> is.null() ) return(NULL)
        
        tibble(
          sample_id = unique(colData(sce_list_per_sample)$sample_id),
          sce = list(sce_list_per_sample)
        )
      },
      pattern = map(sce_list_per_sample),
      iteration = "list"
    ),
    
    # 5. combine to one tibble, separate to different sample
    tar_target(
      sce_list_by_sample_list_split_by_sample,
      bind_rows(sce_list_by_sample_list) |> 
        group_by(sample_id) |>
        tar_group(),
      iteration = "group"
    ),
    
    # 6. cbind sample
    tar_target(
      sample_sce_merged_tbl,
      {
        sce_list_by_sample_list_split_by_sample |>
          summarize(sample_id = dplyr::first(sample_id), sce_list = list(sce), .groups = "drop") |> 
          mutate(
            sce_merged = map(
              sce_list,
              ~ {
                browser()
                xs <- .x     # list of SCE objects
                
                # ---------------------------
                # 1. INTERSECT GENES
                # ---------------------------
                common_genes <- Reduce(intersect, map(xs, rownames))
                xs <- map(xs, ~ .x[common_genes, , drop = FALSE])
                
                # # ---------------------------
                # # 2. UNION colData columns
                # # ---------------------------
                # coldata_union <- Reduce(union, map(xs, ~ colnames(colData(.x))))
                
                xs <- map(xs, function(s) {
                  cd <- colData(s)
                  
                  # # columns missing in this SCE
                  # missing <- setdiff(coldata_union, colnames(cd))
                  # 
                  # # add missing columns as NA
                  # for (m in missing) {
                  #   cd[[m]] <- NA
                  # }
                  # 
                  # # reorder columns
                  # colData(s) <- cd[, coldata_union, drop = FALSE]
                  
                  # REBUILD A SCE ----------------------
                  counts_mat <- assay(s, "counts")
                  counts_mat <- counts_mat[common_genes, , drop = FALSE]
                  
                  SingleCellExperiment::SingleCellExperiment(
                    assays  = list(counts = counts_mat),
                    colData = cd,
                    rowData = S4Vectors::DataFrame(row.names = common_genes)
                  )
           
                })
                
                # ---------------------------
                # 3. CBIND
                # ---------------------------
                do.call(cbind, xs)
              }
            ))
        
      },
      pattern = map(sce_list_by_sample_list_split_by_sample),
      iteration = "list"
    ),
    
    # Save to disk
    tar_target(
      save_h5ad_to_disk,
      save_h5ad(sample_sce_merged_tbl$sample_id, sample_sce_merged_tbl$sce_merged, "/vast/scratch/users/shen.m/htan/hta/09-11-2025/counts/"),
      pattern = map(sample_sce_merged_tbl)
    )
  )
  
}, script = paste0(store, "_target_script.R"), ask = FALSE)


job::job({
  
  tar_make(
    script = paste0(store, "_target_script.R"), 
    store = store, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

tar_workspace(sce_merged_2ed158fc4d9d5a94, store = store, script = paste0(store, "_target_script.R") )
tar_meta(starts_with("sce_list_by_sample_list"), store = store) |> filter(is.na(error)) |> pull(name)
sce_list_by_sample_list = tar_read(sce_list_by_sample_list, store = store) |> bind_rows()

sce_merged = tar_read(sce_merged, store = store) |> bind_rows()
sce_merged |> head() |> pull(sce_merged)




