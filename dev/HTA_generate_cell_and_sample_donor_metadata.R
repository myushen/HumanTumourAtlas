# =============================================================================
# Generate HTAN Cell-Level and Sample-Level Metadata
# =============================================================================
#
# This script generates comprehensive metadata files for Human Tumor Atlas Network
# (HTAN) single-cell RNA sequencing data. It processes h5ad files created by 
# parse_files_to_sce.R and combines them with HTAN metadata files to create unified
# cell-level and sample-level metadata files.
#
# Main Components:
#
# 1. Cell-level metadata extraction:
#    - Extracts cell information from processed h5ad files (from parse_files_to_sce.R)
#    - Extracts cell type information from original synapse h5ad files
#    - Creates key columns:
#      * cell_id: Unique cell identifier (format: barcode___Biospecimen)
#      * sample_id: Biospecimen identifier (or Biospecimen___channel for multi-channel)
#      * file_id_cellNexus_single_cell: File identifier for cellNexus integration
#    - Handles duplicate cells by collapsing cell type annotations
#    - Memory efficient: Only reads colData from h5ad files, not expression data
#
# 2. Sample-level metadata generation:
#    - Processes HTAN files_metadata, samples_metadata, and donors_metadata
#    - Handles multi-channel data (HTAPP: Biospecimen with >4 files get channel suffix)
#    - Handles multi-Biospecimen files (MSK: comma-separated Biospecimen values)
#    - Creates sample_id that matches cell-level metadata
#    - Joins sample and donor-level annotations
#    - Standardizes column names (donor_id, diagnosis_age, center, etc.)
#
# 3. Combined metadata:
#    - Joins cell-level and sample-level metadata
#    - Adds atlas_id for directory hierarchy
#    - Saves as compressed parquet format for efficient storage and querying
#
# Key features:
# - Uses targets package for parallel processing of h5ad files
# - Handles multiple HTAN atlases (HTAPP, BU, MSK) with different data structures
# - Extracts cell type information from synapse h5ad files when available
# - Filters out problematic files (e.g., single_cell_RNAseq_level_4_lung files)
# - Uses DuckDB for efficient querying of large parquet files
# - Includes testing code for cellNexus integration
#
# Output files:
#   - hta_cell_metadata.parquet: Cell-level metadata with cell_id, sample_id, file_id, cell_type
#   - hta_sample_metadata.parquet: Sample-level metadata with all HTAN annotations
#   - hta_metadata.0.1.0.parquet: Combined cell and sample metadata for cellNexus

library(targets)
library(dplyr)
library(duckdb)
# =============================================================================
# Create cell-level metadata
# =============================================================================

store_file_hta_cell_metadata <- "/vast/scratch/users/shen.m/htan/hta_lung_scrnaseq_cell_metadata_target_store"

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
  library(tidySingleCellExperiment)
  library(tidyr)
  
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
  
# Extract cell metadata from each h5ad file
  get_h5ad_cell_metadata <- function(h5ad_file) {
    # Read h5ad file (only colData to be memory efficient)
    sce <- readH5AD(
      file = h5ad_file,
      reader = "R",
      use_hdf5 = TRUE,
      obs = TRUE,  # Read cell metadata
      var = FALSE, # Don't read gene metadata
      raw = FALSE, # Don't read raw data
      layers = FALSE # Don't read layers
    )
    
    # Extract metadata
    col_data <- colData(sce) |>
      as.data.frame() |>
      tibble::rownames_to_column("cell_id") |>
      as_tibble() |> 
      mutate(
        # THERE SHOULD BE ANOTHE COLUMN IN THE METADATA TO SPECIFIY HTA DATA
        file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad")
      ) |>
      dplyr::select(cell_id, sample_id, file_id_cellNexus_single_cell, contains("cell_type"))
    
    col_data
  }
  
  # Extract cell type information from original synapse h5ad metadata
  get_synapse_h5ad_cell_type <- function(h5ad_file) {
    
    map_id <- read.csv("/vast/scratch/users/shen.m/synapse_data/lung/counts/adata_sample_id_htan_id_map.csv")
    
    sce <- readH5AD(
      file = h5ad_file,
      reader = "R",
      use_hdf5 = TRUE,
      obs = TRUE,  # Read cell metadata
      var = FALSE, # Don't read gene metadata
      raw = FALSE, # Don't read raw data
      layers = FALSE # Don't read layers
    )
    
    # THERE ARE SOME H5AD DOESNT HAVE MEANINGFUL METADATA AT ALL
    if (!"patient" %in% colnames(as.data.frame(colData(sce))))
      return(NULL)
    
    wanted <-  grep("patient|cell_type", colnames(colData(sce)), value = TRUE)
    
    colData(sce) <- colData(sce)[ , wanted, drop = FALSE]
    
    sce <- sce |>       tidyr::separate(".cell",
                                 c("sample", "barcode"),
                                 sep = "_(?=[^_]+$)",
                                 remove = FALSE) %>%
      left_join(map_id, by = "sample") %>%
      dplyr::rename(sample_id = sample_HTAN_ID)
    
    col_data <- colData(sce) |>
      as.data.frame() |>
      tibble::rownames_to_column("cell_id") |>
      as_tibble() |> 
      mutate(
        # THERE SHOULD BE ANOTHE COLUMN IN THE METADATA TO SPECIFIY HTA DATA
        file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad")
      ) |> 
      dplyr::select(cell_id, sample_id, file_id_cellNexus_single_cell, contains("cell_type"))

    
    col_data
    
  }
  
  list(
    tar_target(
      h5ad_files,
      # Get all saved h5ad files
      list.files("/vast/scratch/users/shen.m/htan/hta/09-11-2025/counts/", pattern = "\\.h5ad$", full.names = TRUE)
      ),
    
    tar_target(
      synapse_h5ad,
      list.files("/vast/scratch/users/shen.m/synapse_data/lung/counts/", pattern = "\\.h5ad$", full.names = TRUE)
    ),
    tar_target(
      cell_data_list,
      get_h5ad_cell_metadata(h5ad_files),
      pattern = map(h5ad_files)
    ),
    
    tar_target(
      synapse_h5ad_cell_type_df,
      get_synapse_h5ad_cell_type(synapse_h5ad ),
      pattern = map(synapse_h5ad)
    )
  )
}, script = paste0(store_file_hta_cell_metadata, "_target_script.R"), ask = FALSE)


job::job({
  
  tar_make(
    script = paste0(store_file_hta_cell_metadata, "_target_script.R"), 
    store = store_file_hta_cell_metadata, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

# tar_workspace(synapse_h5ad_cell_type_df_1976a09c513942d6, store = store_file_hta_cell_metadata, script = paste0(store_file_hta_cell_metadata, "_target_script.R") )
# debugonce(get_synapse_h5ad_cell_type)
# get_synapse_h5ad_cell_type(synapse_h5ad )

cell_metadata <- tar_read(cell_data_list,  store = store_file_hta_cell_metadata) |> bind_rows()

# Collapse dataframe because cell type data is all over the place
synapse_h5ad_cell_type_df <- tar_read(synapse_h5ad_cell_type_df, store = store_file_hta_cell_metadata ) |> bind_rows() |> 
  filter(!is.na(sample_id)) |> 
  add_count(cell_id, name = "cell_count") %>%
  
  {
    df <- .
    # ---- UNIQUE: keep as-is ----
    unique_cell <- df |> filter(cell_count == 1)
    
    # ---- DUPLICATED: collapse only these ----
    duplicated_cell <- df |> 
      filter(cell_count > 1) |>
      group_by(cell_id, sample_id, file_id_cellNexus_single_cell) |>
      summarize(
        across(everything(), ~ dplyr::first(na.omit(.x))),
        .groups = "drop"
      )
    
    bind_rows(unique_cell, duplicated_cell)
  }

cell_metadata <- cell_metadata |> left_join(synapse_h5ad_cell_type_df |> 
                                              filter(!is.na(sample_id)),
                                            by = c("cell_id", "sample_id", "file_id_cellNexus_single_cell"))


# Save cell metadata
save_directory <- "/vast/scratch/users/shen.m/htan/"

arrow::write_parquet(
  cell_metadata,
  file.path(save_directory, "hta_cell_metadata.parquet")
)

cell_metadata <-  tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_cell_metadata.parquet')")
)



# =============================================================================
# Create Sample, Donor metadata
# =============================================================================
file_metadata <- read.csv("/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv",
           sep = "\t", na.strings = c("NA",""), header = TRUE) |> as_tibble() |> 
  # THIS IS ASSIGNED WRONG TO BIOSPECIMEN. HTAN PHASE1 IS SOOOO COMPLEX!
  filter(Filename != "single_cell_RNAseq_level_4_lung/lung_HTA1_203_332102_ch1_L4.tsv") |>
  
  #filter(File.Format != "hdf5",Atlas.Name != "HTAN MSK") |>
  mutate(
    Filename_basename = basename(Filename),
    full_path = file.path(file_path, Filename_basename)) |> 
  # --- MSK: spread comma-separated Biospecimen into multiple rows ---
  tidyr::separate_rows(
    Biospecimen,
    sep = ", ",
    convert = FALSE
  ) |>
  
  # --- HTAPP: derive channel number (others will just get NA) ---
  mutate(
    channel_number = Filename_basename |>
      str_replace_all("ch(?=[0-9]+)", "channel") |>
      str_extract("channel[0-9]+")
  ) |> # need group size per (Atlas.Name, Biospecimen) to mimic HTAPP logic
  group_by(Atlas.Name, Biospecimen) |>
  mutate(
    n_files_per_biospecimen = dplyr::n(),
    sample_id = case_when(
      # HTAPP: if a Biospecimen has >4 files, append channel; else just Biospecimen
      Atlas.Name == "HTAN HTAPP" & n_files_per_biospecimen > 4 ~
        paste(Biospecimen, channel_number, sep = "___"),
      Atlas.Name == "HTAN HTAPP" ~ Biospecimen,
      
      # BU: sample_id = Biospecimen 
      Atlas.Name == "HTAN BU" ~ Biospecimen,
      
      # MSK: after separate_rows, each Biospecimen is already one per row
      Atlas.Name == "HTAN MSK" ~ Biospecimen,
      
      TRUE ~ Biospecimen
    )
  ) |>
  ungroup() |>
  distinct(sample_id, Biospecimen, Assay, Organ, Atlas.Name, Atlasid) |> 
  mutate(
    file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad")
  ) |> 
  as_tibble()


sample_metadata <-  read.csv(
  "/home/users/allstaff/shen.m/projects/HTAN/samples_metadata_2025_10_21.tsv",
  sep = "\t",
  na.strings = c("NA", ""),
  header = TRUE
)
 
donors_metadata <- read.csv(
  "/home/users/allstaff/shen.m/projects/HTAN/donors_metadata_2025_10_21.tsv",
  sep = "\t",
  na.strings = c("NA", ""),
  header = TRUE
)


# Joining -----------------------------------------------------------------
sample_metadata <- file_metadata |> left_join(
  sample_metadata,
  by = c("Biospecimen" = "HTAN.Biospecimen.ID", "Atlas.Name"),
  copy = TRUE
)

sample_metadata <- sample_metadata |> 
  left_join(
    donors_metadata, 
    by = c("Participant.ID" = "HTAN.Participant.ID", "Atlas.Name", "Publications",
           "Atlas.ID", "Level")
  ) |> 
  dplyr::rename(
    donor_id = Participant.ID,
    diagnosis_age = Age.at.Diagnosis..years.,
    title = Publications,
    center = Atlas.Name,
    center_id = Atlas.ID,
    sex = Gender,
    tissue = Organ,
    url = Protocol.Link,
    primary_diagnosis = Primary.Diagnosis,
    assay = Assay
    
  ) |> 
  mutate(self_reported_ethnicity = ifelse(!is.na(Race) & !is.na(Ethnicity) , paste(Race, Ethnicity, sep = "___"), "unknown"),
         sex = tolower(sex),
         age_days = diagnosis_age * 365,
         tissue = tolower(tissue),
         sex = ifelse(is.na(sex), "unknown", sex),
         assay = ifelse(is.na(assay), "scRNA-seq", assay) # BETTER NOT HARDCODE HERE
         )



# Subset necessary columns ------------------------------------------------
cols <- c(
  "sample_id",
  "donor_id",
  "self_reported_ethnicity",
  "age_days",
  "tissue",
  "sex",
  "assay",
  "file_id_cellNexus_single_cell",
  "atlas_id", 
  #"Biospecimen",
 # "HTAN.Parent.ID",
 # "Source.HTAN.Biospecimen.ID",
  "diagnosis_age",
  "primary_diagnosis",
  "url",
  "center",
  "center_id",
  "title",
  "Last.Known.Disease.Status",
  "Progression.or.Recurrence.Type",
  "Vital.Status",
  "Cause.of.Death",
  "Year.Of.Birth",
  "Year.of.Death",
  "Year.Of.Diagnosis",
  "Tumor.Grade",
  "Progression.or.Recurrence",
  "Tissue.or.Organ.of.Origin"
)

sample_metadata <- sample_metadata |> dplyr::select(any_of(cols))

# sample_metadata |> glimpse()
# sample_metadata |> dplyr::count(is.na(self_reported_ethnicity), is.na(sex), )

# !!!NOTE: sample_metadata here contains those sample_id occur in file_metadata but cant be found/sliced in synapse pre-existing H5AD. 
#       Thus, sample number is not accurate here. This is solved after left join to cell_metadata
sample_metadata |> distinct(sample_id) |> dplyr::count()

# Save Sample metadata
arrow::write_parquet(
  sample_metadata,
  file.path(save_directory, "hta_sample_metadata.parquet")
)


sample_metadata <-  tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_sample_metadata.parquet')")
) 



# Combine cell and sample metadata ----------------------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

cell_metadata <- tbl(
  con,
  sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_cell_metadata.parquet')")
)

sample_metadata <- tbl(
  con,
  sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_sample_metadata.parquet')")
)

cell_sample_metadata <- cell_metadata |> 
  left_join(sample_metadata, by = c("sample_id", "file_id_cellNexus_single_cell"), copy = TRUE) |> 
  mutate(self_reported_ethnicity = ifelse(is.na(self_reported_ethnicity), "unknown", self_reported_ethnicity),
         sex = ifelse(is.na(sex), "unknown", sex),
         assay = ifelse(is.na(assay), "scRNA-seq", assay)) |> # BETTER NOT HARDCODE HERE
  mutate(atlas_id = "hta/09-11-2025") # For directory hierarchy

duckdb_write_parquet <- function(.tbl_sql, path, con) {
  sql_tbl <- dbplyr::sql_render(.tbl_sql)  # render SQL for con
  sql_call <- glue::glue("
    COPY ({sql_tbl}) 
    TO '{path}' 
    (COMPRESSION ZSTD, COMPRESSION_LEVEL 15)
  ")
  dbExecute(con, sql_call)
}

save <- duckdb_write_parquet(
  cell_sample_metadata,
  path = "/vast/scratch/users/shen.m/htan/hta_metadata.0.1.0.parquet",
  con = con
)





# Testing with cellNexus --------------------------------------------------
library(cellNexus)
library(dplyr)
library(testthat)
library(duckdb)
test_that("get_metadata and get_single_cell_experiment return expected SCE for test sample", {
  
  cache <- "/vast/scratch/users/shen.m/htan/"
  save_directory <- tempdir()   # or your defined directory
  ids <- c("HTA1_203_332101.h5ad" ,
           "dc1a2e1504a4b71427b682a6300d02d3___1.h5ad") # One of cellNexus file
  
  cell_metadata <-  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_metadata.0.1.0.parquet')")
  )
  
  # Build test data
  test_data <- cell_metadata |>
    dplyr::filter(file_id_cellNexus_single_cell %in% ids)
  
  # Save parquet test metadata
  test_data |>
    dplyr::collect() |>
    arrow::write_parquet(
      file.path(save_directory, "test_metadata.parquet")
    )
  
  # Read metadata and construct SCE
  sce <- get_metadata(
    cache_directory = cache,
    local_metadata = file.path(save_directory, "test_metadata.parquet")
  ) |>
    dplyr::filter(file_id_cellNexus_single_cell %in% ids ) |>
    cellNexus::get_single_cell_experiment(cache_directory = cache)
  
  # Basic structural checks
  expect_s4_class(sce, "SingleCellExperiment")
  expect_true(ncol(sce) > 0)
  expect_true(nrow(sce) > 0)
  
  # Check assay column exists
  assay_counts <- sce |> dplyr::count(assay)
  expect_true("assay" %in% colnames(colData(sce)))
  expect_s3_class(assay_counts, "data.frame")
  
})

# upload "/vast/scratch/users/shen.m/htan/hta_metadata.0.1.0.parquet" to Nectar cellNexus-metadata/hta_metadata.0.1.0.parquet

