library(dplyr)
library(tibble)
library(glue)
library(purrr)
library(stringr)
library(HPCell)
library(arrow)
library(targets)
library(crew)
library(crew.cluster)
library(duckdb)
library(SummarizedExperiment)

cell_metadata <-  tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_metadata.0.3.0.parquet')")
)
# sce = cell_metadata |>
#   dplyr::filter(sample_id == "HTA1_203_332102"  ) |>
#   cellNexus::get_single_cell_experiment(cache_directory = "/vast/scratch/users/shen.m/htan/")
# 
# sce |> counts() |> as.matrix() |>  hist(ylim=c(0,1e3), breaks = 50)

#breast_cell_metadata <- cell_metadata |> filter(is.na(tissue) | str_detect(tissue, "breast|Breast"))
summary_store = "/vast/scratch/users/shen.m/hta_check_counts_distribution_summary_target_store"
tar_script({
  library(dplyr)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(crew)
  library(crew.cluster)
  library(duckdb)
  
  # Helper to avoid repetition
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
  elastic_500      <- new_elastic("elastic_500",      500, 60 * 24, workers = 3,   crashes_max = 2)
  elastic_300      <- new_elastic("elastic_300",      300, 60 * 24, workers = 8,   crashes_max = 2, backup = elastic_500)
  elastic_160      <- new_elastic("elastic_160",      160, 60 * 24, workers = 8,   crashes_max = 2, backup = elastic_300)
  elastic_120      <- new_elastic("elastic_120",      120, 60 * 4,  workers = 16,  crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
  elastic_80       <- new_elastic("elastic_80",        80, 60 * 4,  workers = 24,  crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
  elastic_40       <- new_elastic("elastic_40",        40, 60 * 4,  workers = 32,  crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
  elastic_20       <- new_elastic("elastic_20",        20, 60 * 4,  workers = 48,  crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
  elastic_10       <- new_elastic("elastic_10",        10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
  elastic_5_minimal <- new_elastic("elastic_5_minimal", 5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)
  
  controllers <- crew_controller_group(
    elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_300, elastic_500, elastic_5_minimal
  )
  
  tar_option_set(
    memory             = "transient",
    garbage_collection = 100,
    storage            = "worker",
    retrieval          = "worker",
    error              = "continue",
    cue                = tar_cue(mode = "thorough"),
    format             = "qs",
    workspace_on_error = TRUE,
    controller         = controllers,
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    )
  )
  
  
  pos_min_med_ratio <- function(x) {
    pos <- x[x > 0]
    min(pos) / median(pos)
  }
  
  pos_min_mean_ratio <- function(x) {
    pos <- x[x > 0]
    min(pos) / mean(pos)
  }
  
  get_positive_mode <- function(x) {
    sort(table(x[x > 0]), decreasing = TRUE)[1] |> names() |> as.numeric()
  }
  
  # Stage 1 – read SCE and return a minimal list with only what downstream needs.
  # Keeping the raw matrix avoids re-reading the file in the metrics stage, while
  # returning a plain list (not a full SCE) keeps the serialised object small.
  read_sce_counts <- function(file) {
    sce        <- readH5AD(file, reader = "R", use_hdf5 = TRUE)
    if (ncol(sce) == 0) return(NULL)
    
    assay_name <- names(sce@assays)[1]
    counts_mat <- assay(sce, assay_name)
    
    list(
      sample_id  = basename(file),
      counts_mat = counts_mat,
      n_cells    = ncol(counts_mat),
      n_genes    = nrow(counts_mat)
    )
  }
  
  # Stage 2 – pure computation on the pre-loaded matrix; no I/O.
  calc_counts_metrics <- function(sce_data) {
    if (is.null(sce_data)) return(NULL)
    
    # Subsample cells to 1e4 if dataset is large
    if (ncol(sce_data$counts_mat) > 5e3) {
      set.seed(42)
      sampled_cols <- sample(ncol(sce_data$counts_mat), 5e3)
      sce_data$counts_mat <- sce_data$counts_mat[, sampled_cols]
    }
    
    counts_vec <- as.numeric(sce_data$counts_mat)
    tol        <- 1e-4
    all_int    <- all(counts_vec == floor(counts_vec), na.rm = TRUE)
    
    tibble::tibble(
      sample_id          = sce_data$sample_id,
      min_val            = min(counts_vec,    na.rm = TRUE),
      median_val         = median(counts_vec, na.rm = TRUE),
      max_val            = max(counts_vec,    na.rm = TRUE),
      counts_gap_min_med = pos_min_med_ratio(counts_vec),
      counts_gap_min_mean= pos_min_mean_ratio(counts_vec),
      positive_mode      = get_positive_mode(counts_vec),
      has_negative       = min(counts_vec, na.rm = TRUE) < 0,
      max_gt_10          = max(counts_vec, na.rm = TRUE) > 10,
      all_integer        = all_int,
      has_floating       = !all_int && all(abs(counts_vec - round(counts_vec)) < tol, na.rm = TRUE),
      n_cells            = sce_data$n_cells,
      n_genes            = sce_data$n_genes
    )
  }
  
  # fx = function(sce) {
  #   counts_vec = assay(sce, "counts") |> as.numeric()
  #   tol        <- 1e-4
  #   all_int    <- all(counts_vec == floor(counts_vec), na.rm = TRUE)
  #   tibble(
  #     sample_id = unique(sce$sample_id),
  #     all_int = all_int
  #   )
  # }
  
  
  list(
    tar_target(
      files,
      list.files(
        "/vast/scratch/users/shen.m/htan/hta_2025/0.3.0/parsed_counts/",
        full.names = TRUE, pattern = "\\.h5ad$"
      ) |>
        # Temporary fix - 0 Biospecimens should not be parsed
        (\(x) x[basename(x) != "0 Biospecimens.h5ad"])(),
      deployment = "main"
    ),
    
    # Stage 1 – I/O-bound; needs memory for the full matrix
    tar_target(
      sce_counts,
      read_sce_counts(files),
      pattern   = map(files),
      iteration = "list",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_10")
      )
    ),
    
    # Stage 2 – CPU-bound, no disk I/O; can run on smaller workers
    tar_target(
      sample_summary_df,
      calc_counts_metrics(sce_counts),
      pattern   = map(sce_counts),
      iteration = "list",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      )
    )
  )
  
}, ask = FALSE, script = glue("{summary_store}/_targets.R"))

job::job({
  
  tar_make(
    reporter = "summary",
    script = glue("{summary_store}/_targets.R"),
    store = glue("{summary_store}/_targets")
  )
  
})

# tar_workspace(sample_summary_df_e8c275e236b86872,script = glue("{summary_store}/_targets.R"),
#               store = glue("{summary_store}/_targets"))
# debugonce(calc_counts_metrics)
# calc_counts_metrics(sce_counts)

sample_summary_df = tar_read(sample_summary_df, store = glue("{summary_store}/_targets")) |>
  dplyr::bind_rows() |>  mutate(max_gt_20 = ifelse(max_val > 20, TRUE, FALSE)) 
# |>
#   mutate(sample_id = stringr::str_remove(sample_id, ".h5ad"))

impute_x_approximate_distribution <- function(df,
                                              counts_gap_threshold,
                                              pos_mode_threshold) {
  df |>
    dplyr::mutate(
      inferred_distribution = dplyr::case_when(
        
        # 0) When counts gap between 0 and next min value >= threshold
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          (counts_gap_min_mean >= counts_gap_threshold) & (positive_mode > pos_mode_threshold) ~ "double_log1p",
        
        # 1) Small counts gap
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          !(
            (counts_gap_min_mean >= counts_gap_threshold) &
              (positive_mode > pos_mode_threshold)
            
          ) ~ "log1p",
        
        # 2) No negatives, has large values
        !has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 3) Large values, integer counts
        !has_negative & max_gt_20 & all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 4) Has negatives, compressed range
        has_negative & !max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 5) Has negatives and large values
        has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 6) No negatives, and all integers
        !has_negative & all_integer ~  "raw_limit_max_to_10",
         
        # fallback
        TRUE ~ NA_character_
      )
    )
}

sample_summary_df = sample_summary_df |> impute_x_approximate_distribution(0.25, 1) |> 
  mutate(count_upper_bound = case_when(
    # 0) When counts gap between 0 and next min value >= 0.25, double log. Max value before exp is 10.
    inferred_distribution == "double_log1p" ~ 10,
    
    # 1) make 10 as max before exp
    inferred_distribution == "log1p" ~ 10,
    
    # 2,3,5) transform algo picks up negative value. should always scale max to 10
    # 4) Has negatives, no large values, no integer, no floating. Counts peak at 10
    inferred_distribution == "raw_limit_max_to_10" ~ 10
    
  )) |>
  # Inverse distribution
  mutate(method_to_apply = case_when(inferred_distribution == "double_log1p" ~ "safe_expm1",
                                     inferred_distribution == "log1p" ~ "expm1",
                                     inferred_distribution == "raw_limit_max_to_10" ~ "identity_with_max_limit"))


sample_summary_df <- sample_summary_df |>
  mutate(file_name = file.path("/vast/scratch/users/shen.m/htan/hta_2025/0.3.0/parsed_counts/", 
                               sample_id),
         feature_thresh = if_else(n_genes > 1e3, 200, floor(500/2e4)*n_genes)) |>
  select(file_name, sample_id, method_to_apply, count_upper_bound, feature_thresh)

sample_summary_df |> arrow::write_parquet("/vast/projects/cellxgene_curated/hta/all_center_sample_summary_df_for_hpcell.parquet", compression = "gzip")
sample_summary_df <- arrow::read_parquet("/vast/projects/cellxgene_curated/hta/all_center_sample_summary_df_for_hpcell.parquet")

sample_names <-
  sample_summary_df |> 
  pull(file_name) |> 
  str_replace("/home/users/allstaff/shen.m/scratch", "/vast/scratch/users/shen.m") |> 
  set_names(sample_summary_df |> pull(sample_id))
functions = sample_summary_df |> pull(method_to_apply)
feature_thresh = sample_summary_df |> pull(feature_thresh)
count_upper_bound = sample_summary_df |> pull(count_upper_bound)

my_store = "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store"

new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 1, backup = NULL) {
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
elastic_300 <- new_elastic("elastic_300", 300, 60 * 24, workers = 4,  crashes_max = 2)
elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2, backup = elastic_300)
elastic_120  <- new_elastic("elastic_120",  120,  60 * 8,  workers = 16, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
elastic_80  <- new_elastic("elastic_80",   80,  60 * 8,  workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)

elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)

# Group for targets (small → large)
controllers <- crew::crew_controller_group(
  elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_300, elastic_5_minimal
)

job::job({
  
  library(HPCell)
  
  sample_names |>
    initialise_hpc(
      store = my_store,
      gene_nomenclature = "ensembl",
      data_container_type = "anndata",
      computing_resources = list(
        elastic_5_minimal, elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_300
      ),
      default_controller = "elastic_10", 
      verbosity = "summary",
      update = "never", 
      #update = "thorough", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    transform_assay(fx = functions, target_output = "sce_transformed", scale_max = count_upper_bound) |>
    
    # # Remove empty outliers based on RNA count threshold per cell
    remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = feature_thresh) |>
    
    # Annotation
    annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |>
    
    # Cell type harmonisation
    celltype_consensus_constructor(target_input = "sce_transformed",
                                   target_output = "cell_type_concensus_tbl", 
                                   annotation_unified_names = c("azimuth", "blueprint", "monaco")) |>
    
    # Alive identification
    remove_dead_scuttle(target_input = "sce_transformed", target_annotation = "cell_type_concensus_tbl",
                        group_by = "cell_type_unified_ensemble") |>
    
    # Doublets identification
    remove_doublets_scDblFinder(target_input = "sce_transformed") |>
    
    # # SCT
    # normalise_abundance_seurat_SCT(target_input = "sce_transformed", factors_to_regress = c(
    #   "subsets_Mito_percent",
    #   "subsets_Ribo_percent")) |>
    # 
    # # Pseudobulk
    # calculate_pseudobulk(target_input = "sce_transformed",
    #                      group_by = "cell_type_unified_ensemble") |>
    
    # # metacell
    # cluster_metacell(target_input = "sce_transformed",  group_by = "cell_type_unified_ensemble") |>
    # 
    # # Cell Chat
    # ligand_receptor_cellchat(target_input = "sce_transformed",
    #                          group_by = "cell_type_unified_ensemble") |>
    
    print()
  
  
})

# Debug
# tar_workspace(cell_type_concensus_tbl_399b737415cfdd3b, store = my_store)
# debugonce(cell_type_ensembl_harmonised)
# cell_type_ensembl_harmonised(sce_transformed, annotation_tbl, NULL, NULL, available_maps = c("azimuth", "blueprint", "monaco"))
# annotation_label_transfer(sce_transformed, empty_tbl, "pbmcref", feature_nomenclature = gene_nomenclature)

tar_progress_branches(store = my_store) |> mutate(pending = branches - skipped - completed)

#' Pipeline for Lightening Annotations in High-Performance Computing Environment
#' 
#' This pipeline is designed to read, process, and "lighten" large annotation tables in an HPC environment.
#' It uses the `targets` package for reproducibility and `crew` for efficient job scheduling on a Slurm cluster.
#' The `lighten_annotation` function selects and processes specific columns from large tables to reduce memory usage.
#' 
#' The pipeline consists of:
#' - **Crew Controllers**: Four tiers of Slurm controllers with varying memory allocations to optimize resource usage.
#' - **Targets**:
#'   - `my_store`: Defines the path to the target storage directory, ensuring all targets use the correct storage location.
#'   - `target_name`: Retrieves metadata to identify branch targets for annotation.
#'   - `annotation_tbl_light`: Applies `lighten_annotation` to process each target name, optimally running with `tier_1` resources.
#' 
#' @libraries:
#'   - `dplyr`, `magrittr`, `tibble`, `targets`, `tarchetypes` for data manipulation and pipeline structure.
#'   - `crew`, `crew.cluster` for parallel computation and cluster scheduling in an HPC environment.
#' 
#' @options:
#'   - Memory settings, garbage collection frequency, and error handling are set to handle large data efficiently.
#'   - The `cue` option is set to `never` for forced target updates if needed.
#'   - `controller` is a group of Slurm controllers to manage computation across memory tiers.
#' 
#' @function `lighten_annotation`: Processes each annotation table target, unnesting and selecting specific columns to reduce data size.
#'
#' @example Usage:
#'   The pipeline script is saved as `/vast/scratch/users/shen.m/lighten_annotation_tbl_target.R` by tar_script and can be run using `tar_make()`.
conflicted::conflict_prefer("filter", "dplyr")
tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
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
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    ),
    controller = controllers, 
    trust_object_timestamps = TRUE
  )
  
  lighten_annotation = function(target_name, my_store ){
    annotation_tbl = tar_read_raw( target_name,  store = my_store )
    if(annotation_tbl |> is.null()) { 
      warning("this annotation is null -> ", target_name)
      return(NULL) 
    }
    
    annotation_tbl |> 
      unnest(blueprint_scores_fine) |> 
      select(.cell, starts_with("sample"), blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), monaco_scores_fine, contains("macro"), contains("CD4") ) |> 
      unnest(monaco_scores_fine) |> 
      select(.cell, starts_with("sample"), blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), contains("macro") , contains("CD4"), contains("helper"), contains("Th")) |> 
      rename(cell_ = .cell)
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store", deployment = "main"),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("annotation_tbl_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    )    ,
    
    tar_target(
      annotation_tbl_light,
      lighten_annotation(target_name, my_store),
      packages = c("dplyr", "tidyr"),
      pattern = map(target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      )
    )
  )
  
  
}, script = "/vast/scratch/users/shen.m/hta_lighten_annotation_tbl_target.R", ask = FALSE)

job::job({
  
  tar_make(
    script = "/vast/scratch/users/shen.m/hta_lighten_annotation_tbl_target.R",
    store = "/vast/scratch/users/shen.m/hta_lighten_annotation_tbl_target", 
    reporter = "summary"
  )
  
})

# Sample metadata
library(arrow)
library(dplyr)
library(duckdb)
library(targets)

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

# Write annotation light
cell_metadata <- 
  tbl(
    con,
    sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_metadata.0.3.0.parquet')")
  )

cell_annotation = 
  tar_read(annotation_tbl_light, store = "/vast/scratch/users/shen.m/hta_lighten_annotation_tbl_target") |> 
  dplyr::rename(
    blueprint_first_labels_fine = blueprint_first.labels.fine, 
    monaco_first_labels_fine = monaco_first.labels.fine, 
    azimuth_predicted_celltype_l2 = azimuth_predicted.celltype.l2
  ) 

cell_annotation = cell_annotation |> mutate(
  blueprint_first_labels_fine = ifelse(is.na(blueprint_first_labels_fine), "Other", blueprint_first_labels_fine),
  monaco_first_labels_fine = ifelse(is.na(monaco_first_labels_fine), "Other", monaco_first_labels_fine),
  azimuth_predicted_celltype_l2=ifelse(is.na(azimuth_predicted_celltype_l2), "Other", azimuth_predicted_celltype_l2))

# # cell_annotation |> arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/annotation_tbl_light.parquet",
# #                                         compression = "zstd")
# cell_annotation |> arrow::write_parquet("~/scratch/cache_temp/annotation_tbl_light.parquet",
#                                         compression = "zstd")

empty_droplet = 
  tar_read(empty_tbl, store = "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store") |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

alive_cells = 
  tar_read(alive_tbl, store = "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store") |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell) |>
  select(-cell_type_unified_ensemble)

doublet_cells =
  tar_read(doublet_tbl, store = "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store") |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

# metacell = 
#   tar_read(metacell_tbl, store = "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store") |> 
#   bind_rows() |> 
#   dplyr::rename(cell_ = cell) |> 
#   dplyr::rename_with(
#     ~ stringr::str_replace(.x, "^gamma", "metacell_"),
#     starts_with("gamma")
#   )

# Save cell type concensus tbl from HPCell output to disk
cell_type_concensus_tbl = tar_read(cell_type_concensus_tbl, store = "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store") |>  
  bind_rows() |> 
  dplyr::rename(cell_ = .cell)

cell_type_concensus_tbl = cell_type_concensus_tbl |> mutate(cell_type_unified_ensemble = 
                                                              ifelse(is.na(cell_type_unified_ensemble),
                                                                     "Unknown",
                                                                     cell_type_unified_ensemble))

cell_metadata_joined = cell_metadata |> 
  left_join(empty_droplet, by = c("cell_id" = "cell_", "sample_id"), copy=TRUE) |>  
  left_join(cell_type_concensus_tbl, by = c("cell_id" = "cell_", "sample_id"),copy=TRUE) |>
  #left_join(cell_annotation, copy=TRUE) |>  
  left_join(alive_cells, by = c("cell_id" = "cell_", "sample_id"),copy=TRUE) |> 
  left_join(doublet_cells, by = c("cell_id" = "cell_", "sample_id"),copy=TRUE)
# |>
#   left_join(metacell, copy=TRUE)

cell_metadata_joined |> filter(is.na(blueprint_first_labels_fine))

# Keep everything as a lazy DuckDB query — no as_tibble(), no in-RAM materialisation.
# coalesce() translates directly to SQL COALESCE
cell_metadata_joined2 = cell_metadata_joined |>
  mutate(
    cell_type_unified_ensemble    = coalesce(cell_type_unified_ensemble,    "Unknown"),
    data_driven_ensemble          = coalesce(data_driven_ensemble,          "Unknown"),
    blueprint_first_labels_fine   = coalesce(blueprint_first_labels_fine,   "Other"),
    monaco_first_labels_fine      = coalesce(monaco_first_labels_fine,      "Other"),
    azimuth_predicted_celltype_l2 = coalesce(azimuth_predicted_celltype_l2, "Other"),
    azimuth                       = coalesce(azimuth,                       "Other"),
    blueprint                     = coalesce(blueprint,                     "Other"),
    monaco                        = coalesce(monaco,                        "Other")
  ) |>
  left_join(
    sample_summary_df |>
      select(sample_id, method_to_apply) |>
      mutate(sample_id = stringr::str_remove(sample_id, ".h5ad")),
    by = "sample_id",
    copy = TRUE
  ) |>
  dplyr::rename(inverse_transform = method_to_apply)

# Unify cell metadata annotations
raw_cols <- cell_metadata_joined2 |> colnames()
pattern_drop <- c(
  grep("^scores", raw_cols, value = TRUE),
  grep("coarse$", raw_cols, value = TRUE)
)
explicit_drop <- c("azimuth",
                   "blueprint",
                   "monaco",
                   "subsets_Mito_sum",
                   "subsets_Mito_detected",
                   "ensemble_joinid",
                   "cell_type_unified",
                   "data_driven_ensemble")
conflicted::conflicts_prefer(base::intersect)
drop_cols <- intersect(unique(c(explicit_drop, pattern_drop)), raw_cols)

# Append pseudobulk count — pb_df stays lazy (same DuckDB connection, no copy needed)
conflicted::conflicts_prefer(dplyr::filter)
pb_df = cell_metadata_joined2 |> 
  filter(!empty_droplet, alive, scDblFinder.class != "doublet") |>
  count(sample_id, cell_type_unified_ensemble, name = ".aggregated_cells")

cell_metadata_joined2 <- cell_metadata_joined2 |> 
  left_join(pb_df, by = c("sample_id", "cell_type_unified_ensemble"))

cell_metadata_joined2 |>
  select(!all_of(drop_cols)) |>
  glimpse()

# Write via DuckDB COPY TO — streams directly to parquet without materialising in R RAM
output_path <- "/vast/projects/cellxgene_curated/hta/metadata_hta.v0.1.0.parquet"
final_sql <- dbplyr::sql_render(cell_metadata_joined2 |> select(!all_of(drop_cols)))
DBI::dbExecute(con, sprintf(
  "COPY (%s) TO '%s' (FORMAT PARQUET, COMPRESSION 'zstd')",
  final_sql, output_path
))

rm(alive_cells, cell_annotation, empty_droplet, cell_type_concensus_tbl, doublet_cells)
gc()

# # Cellchat output
# ligand_receptor_tbl = tar_read(ligand_receptor_tbl, store = "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store") |> bind_rows()

# tar_meta(store = my_store, starts_with("sct")) |> pull(name) |> _[[1]] |>
#   tar_read_raw(store = my_store)


# Save counts, cpm, rank, sct
library(targets)
library(tidyverse)
store_file_cellNexus = "/vast/scratch/users/shen.m/htan/hta/targets_run_normalised_counts_v0.3.0"

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
          distinct(sample_id, file_id_cellNexus_single_cell, file_id_cellNexus_pseudobulk) |> 
          as_tibble(), 
        copy=T
      )
  }
  
  read_target = function(target_name_grouped_by_dataset_id, my_store, target_col, output_col) {
    
    target_col_sym = rlang::sym(target_col)
    output_col_sym = rlang::sym(output_col)
    
    if (any(is.na(target_name_grouped_by_dataset_id[[target_col]]))) {
      return(NULL)
    }
    
    target_name_grouped_by_dataset_id |>
      mutate(
        !!output_col_sym := purrr::map(!!target_col_sym, \(x) {
          tar_read_raw(x, store = my_store)
        })
      )
  }
  
  save_anndata = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    .x = dataset_id_sce |> pull(sce) |> _[[1]]
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    
    .x |> assays() |> names() = "counts"
    
    # Save the experiment data to the specified counts cache directory
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    return(TRUE)  # Indicate successful saving
    
    
  }
  
  save_anndata_cpm = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    # # Parallelise
    dataset_id_sce |> 
      purrr::transpose() |> 
      lapply(
        FUN = function(x) {
          
          .x = x[[ncol(dataset_id_sce)]]
          .y = x[[1]] |> str_remove("\\.h5ad")
          
          # Check if the 'sce' has only one cell (column)
          if(ncol(assay(.x)) == 1) {
            
            # Duplicate the assay to prevent saving errors due to single-column matrices
            my_assay = cbind(assay(.x), assay(.x))
            # Rename the second column to distinguish it
            colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
            
            cd = colData(.x)
            cd = cd |> rbind(cd)
            rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
            
            
            
            .x =  SingleCellExperiment(assay = list( my_assay ) |> set_names(names(assays(.x))[1]), colData = cd) 
          } 
          
          
          # # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
          # .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1], subset_up_to_number_of_cells = 100)
          
          # CALCULATE CPM
          .x =  SingleCellExperiment(assay = list( cpm = calculateCPM(.x, assay.type = names(assays(.x))[1])), colData = colData(.x)) 
          
          # Save the experiment data to the specified counts cache directory
          .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
          
          return(TRUE)  # Indicate successful saving
        }
        
      )
    
    return("saved")
    
  }
  
  process_matrix_in_slices <- function(h5_matrix, output_filepath, output_filepath_temp, chunk_size = 1000) {
    # Load the HDF5 matrix
    n_rows <- dim(h5_matrix)[1]
    n_cols <- dim(h5_matrix)[2]
    
    if (file.exists(output_filepath)) {
      file.remove(output_filepath)
      cat("Existing output file removed.\n")
    }
    if (file.exists(output_filepath_temp)) {
      file.remove(output_filepath_temp)
      cat("Existing output file removed.\n")
    }
    
    # Create an empty list to hold the slices
    slice_list <- list()
    
    # Loop through the matrix in chunks
    for (start_col in seq(1, n_cols, by = chunk_size)) {
      end_col <- min(start_col + chunk_size - 1, n_cols)
      cat("Processing columns", start_col, "to", end_col, "\n")
      
      # Extract a slice of the matrix
      matrix_slice <- as.matrix(h5_matrix[, start_col:end_col, drop=FALSE])
      
      # Calculate ranks for the slice
      ranked_slice <- singscore::rankGenes(matrix_slice)  %>% `-` (1) 
      
      # Convert the ranked slice to sparse format
      sparse_ranked_slice <- as(ranked_slice, "CsparseMatrix")
      
      # Write the slice to the output HDF5 file
      HDF5Array::writeHDF5Array(
        sparse_ranked_slice,
        filepath = output_filepath_temp,
        name = paste0("rank_", start_col, "_to_", end_col),
        as.sparse = TRUE,
        H5type = "H5T_STD_I32LE"
      ) 
      
      # Store the slice name for later binding
      slice_list[[length(slice_list) + 1]] <- paste0("rank_", start_col, "_to_", end_col)
    }
    
    
    slice_list |> map(~HDF5Array::HDF5Array(output_filepath_temp, name =.x)) |> do.call(cbind, args=_)
    
  }
  
  save_rank_per_cell = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)
    
    .x = dataset_id_sce |> pull(sce) |> _[[1]]
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    
    # Check if the 'sce' has only one cell (column)
    if(ncol(assay(.x)) == 1) {
      
      # Duplicate the assay to prevent saving errors due to single-column matrices
      my_assay = cbind(assay(.x), assay(.x))
      # Rename the second column to distinguish it
      colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
      
      cd = colData(.x)
      cd = cd |> rbind(cd)
      rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
      
      
      
      .x =  SingleCellExperiment(assay = list( my_assay ) |> set_names(names(assays(.x))[1]), colData = cd) 
    } 
    
    
    # # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
    # .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1], subset_up_to_number_of_cells = 100)
    
    print("start ranking")
    
    # CALCULATE rank
    rank_assay = 
      .x |>
      assay() |> 
      
      # This because some datasets are still > 1M cells
      process_matrix_in_slices(
        paste(c(cache_directory, "/", .y, "_rank_matrix.HDF5Array"), collapse = ""), 
        paste(c(cache_directory, "/", .y, "_rank_matrix_temp.HDF5Array"), collapse = ""), 
        chunk_size = 1000
      )
    
    print("creating SCE")
    
    .x =  SingleCellExperiment(assay = list( rank = rank_assay), colData = colData(.x)) 
    
    print("saving")
    
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    # Delete the temp file
    file.remove(paste(c(cache_directory, "/", .y, "_rank_matrix_temp.HDF5Array"), collapse = ""))
    
    return(TRUE)  # Indicate successful saving
  }
  
  save_anndata_sct = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    if (is.null(dataset_id_sce)) return(NULL)
    
    .x = dataset_id_sce |> pull(sct) |> _[[1]]
    
    # Fix: check is.null BEFORE ncol() to avoid `argument is of length zero`
    if (is.null(.x) || ncol(.x) == 0) return(NULL)
    
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    
    .x |> assays() |> names() = "sct"
    
    # Wrap save with explicit error logging so the real cause is visible. Strange it shouldnt fail, which passed in debug mode
    tryCatch(
      .x |> save_experiment_data(glue("{cache_directory}/{.y}")),
      error = function(e) {
        message(glue::glue("[save_anndata_sct] FAILED for {.y}: {conditionMessage(e)}"))
        stop(e)
      }
    )
    
    return(TRUE)
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/hta_all_centers_run_hpcell_target_store", deployment = "main"), # MODIFY HERE: HPCell targets store to read SCEs from
    tar_target(cache_directory, "/vast/scratch/users/shen.m/hta_2025/0.3.0", deployment = "main"), # MODIFY HERE: output cache directory for saved anndata files
    tar_target(
      cell_metadata,
      "/vast/projects/cellxgene_curated/hta/metadata_hta.v0.1.0.parquet", # MODIFY HERE: final metadata parquet (should match the COPY TO output above)
      packages = c( "arrow","dplyr","duckdb")
      
    ),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("sce_transformed_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    ),
    
    # tar_target(
    #   sct_target_name,
    #   tar_meta(
    #     starts_with("sct_matrix_"),
    #     store = my_store) |>
    #     filter(type=="branch") |>
    #     pull(name),
    #   deployment = "main"
    # ),
    
    tar_target(
      sample_id_sce_df,
      get_sample_id(target_name, cell_metadata, my_store) |>
        as_tibble() |> dplyr::rename(sce_target_name = target_name),
      packages = "tidySingleCellExperiment",
      pattern = map(target_name)
    ),
    
    # tar_target(
    #   sample_id_sct_df,
    #   get_sample_id(sct_target_name, cell_metadata, my_store) |>
    #     as_tibble() |> dplyr::rename(sct_target_name = target_name),
    #   packages = "tidySingleCellExperiment",
    #   pattern = map(sct_target_name)
    # ),
    
    # join
    tar_target(
      sample_id_target_names_df,
      sample_id_sce_df
      # Uncomment when sct is done pre-calculation
      # |> 
      #   left_join(sample_id_sct_df, by = c("sample_id"), copy=T)
    ),
    
    tar_target(
      target_name_grouped_by_sample_id,
      create_chunks_for_reading_and_saving(sample_id_target_names_df, cell_metadata) |> 
        # 
        # # FOR TESTING PURPOSE ONLY
        # filter(file_id_cellNexus_single_cell %in% c("HTA8_2016_1.h5ad",
        #                                             "HTA1_274_4891101___channel1.h5ad")) |>
        
        group_by(sample_id, file_id_cellNexus_single_cell) |>
        tar_group(),
      iteration = "group",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      ), 
      packages = c("arrow", "duckdb", "dplyr", "glue", "targets")
      
    ),
    
    tar_target(
      sample_id_sce,
      read_target(target_name_grouped_by_sample_id, my_store, "sce_target_name", "sce"),
      pattern = map(target_name_grouped_by_sample_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "glue", "tidyverse", "HPCell", "digest", "scater", "dplyr", "duckdb")
    ),
    
    tar_target(
      saved_anndata,
      save_anndata(sample_id_sce, paste0(cache_directory, "/counts")),
      pattern = map(sample_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "glue", "tidyverse", "HPCell", "digest", "scater", "dplyr", "duckdb")
    )
    
    # 
    # tar_target(
    #   sample_id_sct,
    #   read_target(target_name_grouped_by_sample_id, my_store, "sct_target_name", "sct"),
    #   pattern = map(target_name_grouped_by_sample_id),
    #   packages = c("tidySingleCellExperiment", "SingleCellExperiment", "glue", "tidyverse", "HPCell", "digest", "scater", "dplyr", "duckdb")
    # ),
    # 
    # tar_target(
    #   saved_sample_cpm,
    #   save_anndata_cpm(sample_id_sce, paste0(cache_directory, "/cpm")),
    #   pattern = map(sample_id_sce),
    #   packages = c("tidySingleCellExperiment", "SingleCellExperiment", "glue", "tidyverse", "HPCell", "digest", "scater", "dplyr", "duckdb")
    # ),
    # tar_target(
    #   saved_dataset_rank,
    #   save_rank_per_cell(sample_id_sce, paste0(cache_directory, "/rank")),
    #   pattern = map(sample_id_sce),
    #   packages = c("tidySingleCellExperiment", "SingleCellExperiment", "glue", "tidyverse", "HPCell", "digest", "scater", "dplyr", "duckdb"),
    #   resources = tar_resources(
    #     crew = tar_resources_crew(controller = "elastic_40")
    #   ), 
    # ),
    # tar_target(
    #   saved_sct,
    #   save_anndata_sct(sample_id_sct, paste0(cache_directory, "/sct")),
    #   pattern = map(sample_id_sct),
    #   packages = c("tidySingleCellExperiment", "SingleCellExperiment", "glue", "tidyverse", "HPCell", "digest", "scater", "dplyr", "duckdb")
    # )
    
  )
}, script = paste0(store_file_cellNexus, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(store_file_cellNexus, "_target_script.R"), 
    store = store_file_cellNexus, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

# tar_workspace(saved_pseudobulk_7370e62627ea93e7,store = store_file_cellNexus)
# sample_id_sct_df
# debugonce(save_pseudobulk_anndata)
# save_pseudobulk_anndata(sample_id_pseudobulk, paste0(cache_directory, "/pseudobulk/counts"))
# 
# debugonce(save_anndata_sct)
# save_anndata_sct(sample_id_sct, paste0(cache_directory, "/sct"))
# tar_invalidate(target_name_grouped_by_sample_id, store =store_file_cellNexus)
# 

tar_progress_branches(store = store_file_cellNexus) |> mutate(pending = branches - skipped - completed)

# Then run: ~/git_control/HumanTumourAtlas/dev/prepare_pseudobulk_local_counts_cache.R

# Testing 
library(cellNexus)
library(tidyr)
library(ggplot2)
x = get_metadata(cloud_metadata = NULL, local_metadata = "/vast/projects/cellxgene_curated/hta/metadata_hta.v0.1.0.parquet")
x |> dplyr::count()

qc_summary <- x |>
  mutate(
    category = case_when(
      empty_droplet ~ "Empty droplet",
      !empty_droplet & !alive ~ "Dead cells",
      !empty_droplet & scDblFinder.class == "doublet" ~ "Doublet",
      TRUE ~ "Good cells"
    )
  ) |>
  count(category) |>
  collect() |>
  mutate(pct = n / sum(n) * 100)

# Plot
ggplot(qc_summary, aes(x = reorder(category, -pct), y = pct, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", pct)), vjust = -0.4, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "QC breakdown of HTA",
    x = NULL,
    y = "% of total droplets"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


# # Visualise per tissue QC
# per_tissue_summary <- counts |>
#   filter(
#     (annotation == "empty_droplet" & value == "true")  |
#       (annotation == "alive" & value == "false")          |
#       (annotation == "scDblFinder.class" & value == "doublet")
#   )

sce = x |> keep_quality_cells() |>  # FOR TESTING PURPOSE ONLY
  filter(sample_id %in% c("HTA1_254_571101",
                          "HTA1_141_119201")) |>
  get_single_cell_experiment(assays = c("counts"), cache_directory = "/vast/scratch/users/shen.m/", repository = NULL)


# pb = x |> keep_quality_cells() |> 
#   # FOR TESTING PURPOSE ONLY
#   filter(file_id_cellNexus_single_cell %in% c("HTA1_254_571101.h5ad",
#                                               "HTA1_231_6758320___channel2.h5ad")) |>
#   get_pseudobulk(cache_directory = "/vast/scratch/users/shen.m/htan/", repository = NULL)


# Tissue groups

get_tissue_grouped <- function() {
  list(
    "breast nos" = "Breast",
    "colon nos" = "Colorectal",
    "lung nos" = "Lung",
    "abdomen nos" = "Other and Ill-defined Sites",
    "upper-outer quadrant of breast" = "Breast",
    "blood" = "Bone Marrow",
    "not reported" = "Not Reported",
    "adrenal gland nos" = "Adrenal Gland",
    "ovary" = "Ovary",
    "pancreas nos" = "Pancreas",
    "bone marrow" = "Bone Marrow",
    "skin of trunk, vulva nos, skin of upper limb and shoulder" = "Skin",
    "upper-inner quadrant of breast" = "Breast",
    "lower lobe lung" = "Lung",
    "brain nos" = "Brain",
    "lymph node nos" = "Lymph Nodes",
    "upper lobe lung" = "Lung",
    "unknown" = "Not Reported",
    "liver" = "Liver",
    "cervix uteri" = "Cervix",
    "skin of trunk, unknown, skin of upper limb and shoulder" = "Skin",
    "endometrium" = "Uterus",
    "skin of upper limb and shoulder, skin of trunk, lymph nodes of axilla or arm" = "Skin",
    "lower-outer quadrant of breast" = "Breast",
    "overlapping lesion of breast" = "Breast",
    "mediastinum nos" = "Other and Ill-defined Sites",
    "lower-inner quadrant of breast" = "Breast",
    "tonsil nos" = "Head and Neck",
    "descending colon" = "Colorectal",
    "kidney nos" = "Kidney",
    "ascending colon" = "Colorectal",
    "middle lobe lung" = "Lung",
    "unknown primary site" = "Not Reported",
    "rectum nos" = "Colorectal",
    "skin of upper limb and shoulder" = "Skin",
    "base of tongue nos" = "Head and Neck",
    "retroperitoneum" = "Other and Ill-defined Sites",
    "lower limb nos" = "Other and Ill-defined Sites",
    "long bones of lower limb and associated joints" = "Bone",
    "connective subcutaneous and other soft tissues of thorax" = "Soft Tissue",
    "floor of mouth nos" = "Head and Neck",
    "upper limb nos" = "Other and Ill-defined Sites",
    "rectosigmoid junction" = "Colorectal",
    "paraspinal" = "Other and Ill-defined Sites",
    "cheek mucosa" = "Head and Neck",
    "skin of trunk" = "Skin",
    "skin of lower limb and hip" = "Skin",
    "long bones of upper limb scapula and associated joints" = "Bone",
    "anterior 2/3 of tongue nos" = "Head and Neck",
    "connective subcutaneous and other soft tissues of lower limb and hip" = "Soft Tissue",
    "transverse colon" = "Colorectal",
    "retromolar area" = "Head and Neck",
    "connective subcutaneous and other soft tissues of trunk nos" = "Soft Tissue",
    "larynx nos" = "Head and Neck",
    "sigmoid colon" = "Colorectal",
    "vulva nos" = "Other and Ill-defined Sites",
    "connective subcutaneous and other soft tissues nos" = "Soft Tissue",
    "connective subcutaneous and other soft tissues of abdomen" = "Soft Tissue",
    "cerebellum nos" = "Brain",
    "parietal lobe" = "Brain",
    "temporal lobe" = "Brain"
  ) |>
    enframe(name = "tissue", value = "tissue_groups") |>
    unnest(tissue_groups) |>
    distinct()
}

x |> left_join(get_tissue_lookup(), copy=T) |> count(tissue_groups) |> arrange(desc(n))


