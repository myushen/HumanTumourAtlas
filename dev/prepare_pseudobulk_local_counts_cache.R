# Save cpm, rank, sct
library(targets)
library(tidyverse)
store_file_cellNexus = "/vast/scratch/users/shen.m/htan/hta/targets_run_pseudobulk"

tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  library(glue)
  
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
    cue = tar_cue(mode = "never"),
    format = "qs",
    #debug = "dataset_id_sct_ea377f6e2d0ae2b7",
    workspace_on_error = TRUE,
    controller = controllers, 
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    ) 
  )
  
  get_sample_id = function(target_name, file_id_db_file, my_store){
    # Try reading the target safely (for some failing targets)
    sce = tryCatch(
      tar_read_raw(target_name, store = my_store),
      error = function(e) return(NULL)
    )
    
    # Still need to catch target_name
    if(sce |> is.null()) return(tibble(sample_id = NA_character_, 
                                       target_name= !!target_name))
    
    sce |> 
      
      distinct(sample_id) |> mutate(target_name = !!target_name)
    
  }
  
  create_chunks_for_reading_and_saving = function(dataset_id_sample_id, cell_metadata){
    
    # Solve sample_id mismatches because some end with .h5ad suffix while others dont 
    dataset_id_sample_id |> 
      
      left_join(
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql(glue("SELECT * FROM read_parquet('{cell_metadata}')"))
        )   |> 
          distinct(sample_id, cell_type_unified_ensemble, file_id_cellNexus_pseudobulk) |> 
          as_tibble(), 
        copy=T
      )
  }
  
  # read_target = function(target_name_grouped_by_dataset_id, my_store, target_col, output_col) {
  #   
  #   target_col_sym = rlang::sym(target_col)
  #   output_col_sym = rlang::sym(output_col)
  #   
  #   if (any(is.na(target_name_grouped_by_dataset_id[[target_col]]))) {
  #     return(NULL)
  #   }
  #   
  #   target_name_grouped_by_dataset_id |>
  #     mutate(
  #       !!output_col_sym := purrr::map(!!target_col_sym, \(x) {
  #         tar_read_raw(x, store = my_store)
  #       })
  #     )
  # }
  # 
  save_anndata = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    .x = dataset_id_sce |> pull(sce) |> _[[1]]
    .y = dataset_id_sce |> pull(file_id_cellNexus_pseudobulk) |> _[[1]] |> str_remove("\\.h5ad")
    
    .x |> assays() |> names() = "counts"
    
    # Drop list-type columns in colData 
    cd <- colData(.x)
    is_list_col <- vapply(cd, is.list, logical(1))
    colData(.x) <- cd[, !is_list_col, drop = FALSE]
    
    # Check if there is a memory issue 
    assays(.x) <- assays(.x) |> map(DelayedArray::realize)
    
    # Save the experiment data to the specified counts cache directory
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    return(TRUE)  # Indicate successful saving
    
  }
  
  cbind_sce_by_dataset_id = function(target_name_grouped_by_dataset_id, 
                                     file_id_db_file, my_store){
    
    my_cell_type = unique(target_name_grouped_by_dataset_id$cell_type_unified_ensemble)
    
    file_id_db = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{file_id_db_file}')"))
      ) |> 
      dplyr::filter(cell_type_unified_ensemble %in% my_cell_type) |>
      select(sample_id, cell_type_unified_ensemble,
             file_id_cellNexus_pseudobulk) 
    
    
    file_id_db = 
      target_name_grouped_by_dataset_id |> 
      left_join(file_id_db, copy = TRUE)
    
    
    # Parallelise
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)) -1
    # Respect R CMD CHECK core limit if set
    if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
      cores <- min(cores, 2L)
    }
    bp <-   MulticoreParam(workers = cores, progressbar = TRUE)
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    sce_df = 
      file_id_db |> 
      mutate(cell_id = paste(sample_id, cell_type_unified_ensemble, sep = "___")) |>
      nest(cells = cell_id) |> 
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sce = bplapply(
          target_name,
          FUN = function(x) tar_read_raw(x, store = my_store),
          BPPARAM = bp
        )
      ) |>
      
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sce, is.null)) |> 
      
      mutate(sce = map2(sce, cells, ~ .x |>
                          filter(.cell %in% .y$cell_id),
                        
                        .progress = TRUE))
    
    
    
    if(nrow(sce_df) == 0) {
      warning("this chunk has no rows for somereason.")
      return(NULL)
    }
    
    sce_df = sce_df |> 
      
      # THIS SHOULD HAVE BEEN DONE IN THE TRANFORM HPCell
      mutate(sce = map(sce, ~  SingleCellExperiment(assay = assays(.x), colData = colData(.x)) ))
    
    
    # Extra Step 1: Harmonize colData columns - Avoid column name mismatch, force cbind
    all_col_names <- sce_df$sce %>%
      map(~colnames(colData(.x))) %>% 
      unlist() %>% 
      unique()
    
    # Extra Step 2: Standardize colData to have the same columns in each SCE
    sce_df$sce <- map(sce_df$sce, function(sce) {
      current_cols <- colnames(colData(sce))
      missing_cols <- setdiff(all_col_names, current_cols)
      
      if (length(missing_cols) > 0) {
        
        # Fill missing colData columns with NA
        for (col in missing_cols) {
          # Handle sce with empty cells
          if (ncol(sce) == 0)  colData(sce)[, col] <- character(0)
          else if (ncol(sce) > 0) colData(sce)[, col] <- NA
        }
      }
      
      # Ensure the order of columns matches
      colData(sce) <- colData(sce)[, all_col_names]
      return(sce)
    })
    
    sce_df |>
      
      # Step 5: Combine all 'sce' objects within each group into a single 'sce' object
      group_by(file_id_cellNexus_pseudobulk) |>
      summarise( sce =  list(do.call(cbind, args = sce) ) ) 
    
  }
  
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/hta_lung_run_hpcell_target_store", deployment = "main"), # MODIFY HERE: HPCell targets store to read SCEs from
    tar_target(cache_directory, "/vast/scratch/users/shen.m/htan/hta_2025/0.1.0/pseudobulk", deployment = "main"), # MODIFY HERE: output cache directory for saved anndata files
    tar_target(
      cell_metadata,
      "/vast/projects/cellxgene_curated/hta/metadata_hta_lung.v0.1.0.parquet", # MODIFY HERE: final metadata parquet (should match the COPY TO output above)
      packages = c( "arrow","dplyr","duckdb")
      
    ),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("pseudobulk_se_iterated_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    ),
    
    tar_target(
      sample_id_sce_df,
      get_sample_id(target_name, cell_metadata, my_store),
      packages = "tidySingleCellExperiment",
      pattern = map(target_name)
    ),
    
    tar_target(
      target_name_grouped_by_sample_id,
      create_chunks_for_reading_and_saving(sample_id_sce_df, cell_metadata) |> 

        # # FOR TESTING PURPOSE ONLY
        # filter(file_id_cellNexus_pseudobulk %in% c("HTA8_2016_1.h5ad",
        #                                             "HTA1_274_4891101___channel1.h5ad")) |>
        
        group_by(sample_id, file_id_cellNexus_pseudobulk) |>
        tar_group(),
      iteration = "group",
      packages = c("arrow", "duckdb", "dplyr", "glue", "targets")
      
    ),
    
    tar_target(
      sample_id_sce,
      cbind_sce_by_dataset_id(target_name_grouped_by_sample_id, cell_metadata, my_store),
      pattern = map(target_name_grouped_by_sample_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "scater", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly"),
    ),
    
    tar_target(
      saved_pseudobulk,
      save_anndata(sample_id_sce, paste0(cache_directory, "/counts")),
      pattern = map(sample_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "glue", "tidyverse", "HPCell", "digest", "scater", "dplyr", "duckdb")
    )
    
  )
}, script = paste0(store_file_cellNexus, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(store_file_cellNexus, "_target_script.R"), 
    store = store_file_cellNexus, 
    reporter = "summary" #, callr_function = NULL
  )
  
})
# tar_workspace(sample_id_sce_57f89462305fac5b,store=store_file_cellNexus)
# debugonce(cbind_sce_by_dataset_id)
# cbind_sce_by_dataset_id(target_name_grouped_by_sample_id, cell_metadata, my_store)
# tar_invalidate(target_name_grouped_by_sample_id, store =store_file_cellNexus)

x = get_metadata(cloud_metadata = NULL, local_metadata = "/vast/projects/cellxgene_curated/hta/metadata_hta_lung.v0.1.0.parquet")

pb= x |> keep_quality_cells() |> 
#   # FOR TESTING PURPOSE ONLY
#   filter(file_id_cellNexus_pseudobulk %in% c("HTA8_2007_1.h5ad")) |>
# #                                              "HTA1_274_4891101___channel1.h5ad")) |>
  get_pseudobulk(cache_directory = "/vast/scratch/users/shen.m/htan/", repository = NULL)


library(AnnotationDbi)
library(org.Hs.eg.db)
plot_pseudobulk_boxplot = function(pb, top_n = 20, cell_type_col) {
  
  lc = assay(pb, "logcounts")
  
  # Convert Ensembl IDs to gene symbols
  gene_map = AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = rownames(lc),
    column  = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Replace rownames, keep Ensembl ID as fallback if no symbol found
  rownames(lc) = ifelse(is.na(gene_map[rownames(lc)]), rownames(lc), gene_map[rownames(lc)])
  
  # Top variable genes
  top_genes = names(sort(rowVars(lc), decreasing = TRUE)[1:top_n])
  
  # Tidy format
  plot_df = lc[top_genes, ] |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "logcounts") |>
    left_join(
      colData(pb) |> as.data.frame() |> tibble::rownames_to_column("sample") |>
        dplyr::select(sample, cell_type = all_of(cell_type_col)),
      by = "sample"
    )
  
  # Plot
  ggplot(plot_df, aes(x = cell_type, y = logcounts, fill = cell_type)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_fill_manual(
      values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_distinct(plot_df$cell_type))
    ) +
    labs(
      x = "Cell type", y = "logcounts",
      title = glue::glue("Top {top_n} variable genes across cell types")
    ) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(face = "italic")
    )
}

# normalise 
pb = scuttle::logNormCounts(pb) 

# Plot
plot_pseudobulk_boxplot(pb, top_n = 10, cell_type_col = "cell_type_unified_ensemble")

