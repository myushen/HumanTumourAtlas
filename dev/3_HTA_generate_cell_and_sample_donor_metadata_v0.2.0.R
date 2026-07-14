# =============================================================================
# Generate HTAN Cell-Level and Sample-Level Metadata — v0.2.0
# =============================================================================
#
# Reads the h5ad files produced by 2_HTAN_parse_all_centers_targets.R and
# combines them with HTAN metadata files to create unified cell-level and
# sample-level metadata files for cellNexus.
#
# Key differences from v0.1.0:
#   - Source h5ad directory : 0.2.0/counts/  (lung + breast combined)
#   - Cell index store      : breast_all_centers_target_store  (all 7 centers)
#   - New centers           : OHSU, WUSTL, Duke, DFCI  (breast)
#   - DFCI biospecimens use comma-separated values in files_metadata → separate_rows
#   - Cell type columns are read directly from h5ad colData; no separate
#     synapse h5ad extraction step required
#
# Main Components:
#
# 1. Cell-level metadata extraction (via targets for parallelism):
#    - Reads colData only from each h5ad (memory efficient)
#    - Creates cell_id, sample_id, file_id_cellNexus_single_cell/pseudobulk
#    - Preserves any cell_type_* columns present in colData
#
# 2. Cell barcode index map:
#    - Reads cell_index_df targets from breast_all_centers_target_store
#    - Maps numeric cell_id → original barcode (.cell) per center
#
# 3. Sample-level metadata generation:
#    - Processes combined lung+breast files_metadata, samples_metadata, donors_metadata
#    - Handles HTAPP channel splitting, MSK/DFCI comma-separated Biospecimen
#    - Standardises column names (donor_id, diagnosis_age, center, sex, etc.)
#
# 4. Combined metadata:
#    - Joins cell → barcode map → sample/donor annotations
#    - Adds atlas_id = "hta_2026/0.3.0"
#    - Saves as ZSTD-compressed parquet
#
# Output files:
#   - hta_cell_metadata_v020.parquet
#   - hta_sample_metadata_v020.parquet
#   - hta_metadata.0.2.0.parquet

library(targets)
library(dplyr)
library(duckdb)
library(stringr)

# =============================================================================
# 1. Cell-level metadata extraction
# =============================================================================

store_file_hta_cell_metadata <- "/vast/scratch/users/shen.m/htan/hta_v030_cell_metadata_target_store"

tar_script({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(targets)
  library(crew)
  library(crew.cluster)
  
  new_elastic <- function(name, mem_gb, time_min, workers, crashes_max,
                          cpus_per_task = 1, backup = NULL) {
    crew_controller_slurm(
      name            = name,
      workers         = workers,
      crashes_max     = crashes_max,
      seconds_idle    = 30,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = mem_gb,
        cpus_per_task             = cpus_per_task,
        time_minutes              = time_min
      ),
      backup = backup
    )
  }
  
  elastic_160        <- new_elastic("elastic_160",       160, 60 * 24, workers = 8,   crashes_max = 2, cpus_per_task = 2)
  elastic_80         <- new_elastic("elastic_80",         80, 60 * 4,  workers = 24,  crashes_max = 1, backup = elastic_160)
  elastic_40         <- new_elastic("elastic_40",         40, 60 * 4,  workers = 32,  crashes_max = 1, backup = elastic_80)
  elastic_20         <- new_elastic("elastic_20",         20, 60 * 4,  workers = 48,  crashes_max = 1, backup = elastic_40)
  elastic_10         <- new_elastic("elastic_10",         10, 60 * 4,  workers = 150, crashes_max = 2, backup = elastic_20)
  elastic_5_minimal  <- new_elastic("elastic_5_minimal",   5, 60 * 4,  workers = 300, crashes_max = 2, backup = elastic_10)
  
  controllers <- crew_controller_group(
    elastic_5_minimal, elastic_10, elastic_20, elastic_40, elastic_80, elastic_160
  )
  
  tar_option_set(
    memory                  = "transient",
    garbage_collection      = 100,
    storage                 = "worker",
    retrieval               = "worker",
    error                   = "continue",
    cue                     = tar_cue(mode = "thorough"),
    format                  = "qs",
    workspace_on_error      = TRUE,
    controller              = controllers,
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    )
  )
  
  # Extract colData only — memory efficient
  get_h5ad_cell_metadata <- function(h5ad_file) {
    sce <- zellkonverter::readH5AD(
      file     = h5ad_file,
      reader   = "R",
      use_hdf5 = TRUE,
      obs      = TRUE,
      var      = FALSE,
      raw      = FALSE,
      layers   = FALSE
    )
    
    colData(sce) |>
      as.data.frame() |>
      tibble::rownames_to_column("cell_id") |>
      tibble::as_tibble() |>
      mutate(
        file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad"),
        file_id_cellNexus_pseudobulk  = file_id_cellNexus_single_cell
      ) |>
      dplyr::select(
        cell_id, sample_id, contains("barcode"),
        file_id_cellNexus_single_cell, file_id_cellNexus_pseudobulk,
        dplyr::contains("cell_type")
      )
  }
  
  list(
    tar_target(
      h5ad_files,
      list.files(
        "/vast/scratch/users/shen.m/htan/hta_2026/0.3.0/parsed_counts/",
        pattern   = "\\.h5ad$",
        full.names = TRUE
      )
    ),
    
    tar_target(
      cell_data_list,
      get_h5ad_cell_metadata(h5ad_files),
      pattern = map(h5ad_files)
    )
  )
}, script = paste0(store_file_hta_cell_metadata, "_target_script.R"), ask = FALSE)


job::job({
  tar_make(
    script   = paste0(store_file_hta_cell_metadata, "_target_script.R"),
    store    = store_file_hta_cell_metadata,
    reporter = "summary"
  )
})


# =============================================================================
# 2. Read cell metadata and build barcode index map
# =============================================================================

cell_metadata <- tar_read(cell_data_list, store = store_file_hta_cell_metadata) |>
  dplyr::bind_rows() |>
  mutate(cell_id = as.numeric(cell_id)) |>
  mutate(barcode = coalesce(barcode, Barcode)) |>
  select(-Barcode) |>
  # Temporary fix: Why 0 biospecimens were parsed? Need to investigate.
  filter(sample_id != "0 Biospecimens")

# This probably does not matter, because all index were assigned in 2_HTAN_parse_all_centers_targets.R
# Store produced by 2_HTAN_parse_all_centers_targets.R  (organ = "breast",
# but the pipeline covers BOTH lung and breast)
# store_all_centers <- "/vast/scratch/users/shen.m/htan/breast_all_centers_target_store"
# 
# # --- HTAPP ---------------------------------------------------------------
# htapp_index <- tar_read(htapp_cell_index_full_df, store = store_all_centers) |>
#   dplyr::rename(.cell = NAME)
# 
# # --- BU ------------------------------------------------------------------
# bu_index <- tar_read(bu_cell_index_df, store = store_all_centers) |>
#   dplyr::rename(.cell = NAME, sample_id = SampleID)
# 
# # --- MSK -----------------------------------------------------------------
# msk_index <- tar_read(msk_cell_index_df, store = store_all_centers) |>
#   dplyr::select(.cell, sample_id, cell_index)
# 
# # --- OHSU ----------------------------------------------------------------
# ohsu_index <- tar_read(ohsu_cell_index_df, store = store_all_centers) |>
#   dplyr::rename(.cell = NAME)
# 
# # --- WUSTL  (patterned target → list; bind_rows needed) -----------------
# wustl_index <- tar_read(wustl_cell_index_df, store = store_all_centers) |>
#   dplyr::bind_rows() |>
#   dplyr::rename(.cell = NAME)
# 
# # --- Duke ----------------------------------------------------------------
# duke_index <- tar_read(duke_cell_index_df, store = store_all_centers) |>
#   dplyr::rename(.cell = NAME)
# 
# # --- DFCI ----------------------------------------------------------------
# dfci_index <- tar_read(dfci_cell_index_df, store = store_all_centers) |>
#   dplyr::rename(.cell = NAME)
# 
# # Combined index: numeric cell_index → original barcode (.cell) per sample
# cell_index_map <- dplyr::bind_rows(
#   htapp_index,
#   bu_index,
#   msk_index,
#   ohsu_index,
#   wustl_index,
#   duke_index,
#   dfci_index
# ) |>
#   mutate(cell_index = as.numeric(cell_index))
# 
# # Attach original barcodes to cell metadata. This probably do not matter because cell index has been assigned to h5ad
# cell_metadata <- cell_metadata |>
#   dplyr::left_join(cell_index_map, by = c("cell_id" = ".cell", "sample_id"))

cell_metadata |> arrow::write_parquet("/vast/scratch/users/shen.m/htan/cell_metadata_v0.3.0.parquet")

dplyr::tbl(
  DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  dplyr::sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/cell_metadata_v0.3.0.parquet')")
)

rm(list = ls(pattern = "index$"))
gc()

# =============================================================================
# 3. Sample and donor metadata
# =============================================================================

save_directory   <- "/vast/scratch/users/shen.m/htan/"
files_meta_path  <- "/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/files_metadata_scRNAseq_synapse_level3_4.tsv"
samples_meta_path <- "/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/samples_metadata_scRNAseq_synapse_level3_4.tsv"
donors_meta_path  <- "/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/donors_metadata_scRNAseq_synapse_level3_4.tsv"

file_metadata_raw <- read.csv(
  files_meta_path,
  sep        = "\t",
  na.strings = c("NA", ""),
  header     = TRUE
) |>
  tibble::as_tibble() |>
  filter(Filename != "single_cell_RNAseq_level_4_lung/lung_HTA1_203_332102_ch1_L4.tsv") |>
  mutate(
    Filename_basename = basename(Filename),
    # run_dir used by CHOP L3 multi-run logic below
    run_dir           = basename(dirname(Filename)),
    # HTAPP: library_channel mirrors the targets pipeline derivation exactly:
    #   library_id      = segment immediately before _channel{N}
    #   library_channel = paste(library_id, channel_number, sep = "_")
    filename_renamed  = Filename_basename |>
      stringr::str_replace_all("ch(?=[0-9]+)", "channel"),
    channel_number    = filename_renamed |>
      stringr::str_extract("channel[0-9]+"),
    library_id        = filename_renamed |>
      stringr::str_extract("[^_]+(?=_channel[0-9]+)"),
    library_channel   = dplyr::if_else(
      !is.na(library_id) & !is.na(channel_number),
      stringr::str_c(library_id, "_", channel_number),
      NA_character_
    )
  ) |>
  # MSK and DFCI: comma-separated Biospecimen IDs → one row per biospecimen
  tidyr::separate_rows(Biospecimen, sep = ",\\s*", convert = FALSE)

# Pre-compute which CHOP L3 biospecimens span multiple run directories.
# Applies the same file filters as chop_l3_metadata in the targets script.
chop_l3_multirun_bios <- file_metadata_raw |>
  dplyr::filter(
    Atlas.Name == "HTAN CHOP", Level == "Level 3",
    !grepl(",", Biospecimen), Biospecimen != "0 Biospecimens",
    !stringr::str_detect(Filename, "_bak/"),
    !stringr::str_detect(Filename, "snRNA"),
    Filename_basename %in% c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  ) |>
  dplyr::group_by(Biospecimen) |>
  dplyr::summarise(n_runs = dplyr::n_distinct(run_dir), .groups = "drop") |>
  dplyr::filter(n_runs > 1L) |>
  dplyr::pull(Biospecimen)

file_metadata <- file_metadata_raw |>
  dplyr::group_by(Atlas.Name, Biospecimen) |>
  dplyr::mutate(
    n_files_per_biospecimen = dplyr::n(),
    sample_id = dplyr::case_when(
      # HTAPP: biospecimen with >4 files → library_channel suffix
      # Matches: paste(Biospecimen, library_channel, sep="___") in targets
      Atlas.Name == "HTAN HTAPP" & n_files_per_biospecimen > 4 ~
        paste(Biospecimen, library_channel, sep = "___"),
      # CHOP L3 multi-run → run_dir suffix, matching chop_l3_metadata target
      Atlas.Name == "HTAN CHOP" & Level == "Level 3" &
        Biospecimen %in% chop_l3_multirun_bios ~
        paste(Biospecimen, run_dir, sep = "___"),
      # All other centers: sample_id = Biospecimen
      TRUE ~ Biospecimen
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::distinct(sample_id, Biospecimen, Assay, Organ, Atlas.Name, Atlasid) |>
  dplyr::mutate(
    file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad"),
    file_id_cellNexus_pseudobulk  = file_id_cellNexus_single_cell
  )

#  (Investigate further) Biospecimen, center id not parsed
# existing_ids <- list.files("/vast/scratch/users/shen.m/htan/hta_2026/0.3.0/parsed_counts/")
# file_metadata |>
#   filter(file_id_cellNexus_single_cell %in% existing_ids) |>
#   dplyr::count(Atlas.Name)

sample_metadata <- read.csv(
  samples_meta_path,
  sep        = "\t",
  na.strings = c("NA", ""),
  header     = TRUE
) |>
  dplyr::select(where(~ !all(is.na(.))))

donors_metadata <- read.csv(
  donors_meta_path,
  sep        = "\t",
  na.strings = c("NA", ""),
  header     = TRUE
) |>
  dplyr::select(where(~ !all(is.na(.))))

# Join file → sample → donor
sample_metadata <- file_metadata |>
  filter(Biospecimen |> str_detect("HTA")) |>
  dplyr::left_join(
    sample_metadata,
    by = c("Biospecimen" = "HTAN.Biospecimen.ID", "Atlas.Name"),
    copy = TRUE
  ) |>
  dplyr::left_join(
    donors_metadata,
    by = c("Participant.ID" = "HTAN.Participant.ID", "Atlas.Name", "Atlas.ID", "Level")
  ) |>
  dplyr::rename(
    donor_id         = Participant.ID,
    diagnosis_age    = Age.at.Diagnosis..years.,
    center           = Atlas.Name,
    center_id        = Atlas.ID,
    sex              = Gender,
    tissue           = Organ,
    url              = Protocol.Link,
    primary_diagnosis = Primary.Diagnosis,
    assay            = Assay
  ) |>
  mutate(
    self_reported_ethnicity = dplyr::if_else(
      !is.na(Race) & !is.na(Ethnicity),
      paste(Race, Ethnicity, sep = "___"),
      "unknown"
    ),
    sex      = tolower(sex),
    age_days = diagnosis_age * 365,
    tissue   = tolower(tissue),
    sex      = dplyr::if_else(is.na(sex),   "unknown",   sex),
    assay    = dplyr::if_else(is.na(assay), "scRNA-seq", assay)
  ) |>
  dplyr::select(-dplyr::any_of(c("Synapse.ID.x", "Synapse.ID.y"))) |> 
  filter(!sample_id |> str_detect("NA")) |>
  group_by(sample_id) %>%             
  slice_min(nchar(tissue), n = 1) %>% # some center collapse tissues to one row separated by comma. choose the short length. 
  ungroup()

sample_metadata |> dplyr::distinct(sample_id) |> dplyr::count()
sample_metadata|>dplyr::filter(is.na(tissue)) |> dplyr::count(center) # unannotated in the source file metadata

arrow::write_parquet(
  sample_metadata,
  file.path(save_directory, "hta_sample_metadata_v0.3.0.parquet")
)


# =============================================================================
# 4. Combine cell and sample metadata
# =============================================================================

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")

cell_tbl <- dplyr::tbl(
  con,
  dplyr::sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/cell_metadata_v0.3.0.parquet')")
)

sample_tbl <- dplyr::tbl(
  con,
  dplyr::sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_sample_metadata_v0.3.0.parquet')")
)

cell_sample_metadata <- cell_tbl |>
  dplyr::left_join(
    sample_tbl,
    by = c("sample_id", "file_id_cellNexus_single_cell", "file_id_cellNexus_pseudobulk"),
    copy = TRUE
  ) |>
  mutate(
    self_reported_ethnicity = dplyr::if_else(is.na(self_reported_ethnicity), "unknown", self_reported_ethnicity),
    sex   = dplyr::if_else(is.na(sex),   "unknown",   sex),
    assay = dplyr::if_else(is.na(assay), "scRNA-seq", assay),
    atlas_id = "hta_2026/0.3.0"
  )

duckdb_write_parquet <- function(.tbl_sql, path, con) {
  sql_tbl  <- dbplyr::sql_render(.tbl_sql)
  sql_call <- glue::glue("
    COPY ({sql_tbl})
    TO '{path}'
    (COMPRESSION ZSTD, COMPRESSION_LEVEL 15)
  ")
  DBI::dbExecute(con, sql_call)
}

duckdb_write_parquet(
  cell_sample_metadata,
  path = "/vast/scratch/users/shen.m/htan/hta_metadata.0.3.0.parquet",
  con  = con
)

DBI::dbDisconnect(con, shutdown = TRUE)


# =============================================================================
# Testing with cellNexus
# =============================================================================

library(cellNexus)
library(testthat)
library(tidySingleCellExperiment)

test_that("get_metadata and get_single_cell_experiment return expected SCE", {
  
  save_directory <- tempdir()
  
  # Pick one sample from each major organ to test breadth
  ids <- c(
    "HTA1_213_6752601___channel1.h5ad",                      # HTAPP lung
    "dc1a2e1504a4b71427b682a6300d02d3___1.h5ad"  # MSK lung
  )
  
  combined_meta <- dplyr::tbl(
    DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    dplyr::sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_metadata.0.2.0.parquet')")
  )
  
  test_data <- combined_meta |>
    dplyr::filter(file_id_cellNexus_single_cell %in% ids) |>
    dplyr::collect()
  
  arrow::write_parquet(test_data, file.path(save_directory, "test_metadata.parquet"))
  
  sce <- cellNexus::get_metadata(
    local_metadata = file.path(save_directory, "test_metadata.parquet")
  ) |>
    dplyr::filter(file_id_cellNexus_single_cell %in% ids) |>
    cellNexus::get_single_cell_experiment(repository = NULL,
                                          cache_directory = "/vast/scratch/users/shen.m/htan/")
  
  expect_s4_class(sce, "SingleCellExperiment")
  expect_true(ncol(sce) > 0)
  expect_true(nrow(sce) > 0)
  expect_true("assay" %in% colnames(colData(sce)))
})

# upload "/vast/scratch/users/shen.m/htan/hta_metadata.0.2.0.parquet"
#     to Nectar cellNexus-metadata/hta_metadata.0.2.0.parquet
