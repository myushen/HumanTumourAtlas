# =============================================================================
# Unified HTAN targets pipeline — dispatches per center using the registry
#
# Usage:
#   Set `organ` and the three metadata paths below, then source this file.
#   The same script works for Lung, Breast, or any future organ.
#
# Outputs: one .h5ad per HTAN biospecimen ID, saved to `output_dir`.
# =============================================================================

library(targets)
library(dplyr)

# Bind files_metadata
read.csv("/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv",
         sep = "\t", na.strings = c("NA",""), header = TRUE) |> bind_rows(
           read.csv("/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/breast/files_metadata.tsv",
                    sep = "\t", na.strings = c("NA",""), header = TRUE) 
         ) |> as_tibble() |>
  write.table(file = "/home/users/allstaff/shen.m/projects/HTAN/files_metadata_lung_breastcombined.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = "")


# ── Configure these four lines per run ───────────────────────────────────────
organ            <- "breast"          # used only for naming the store
files_meta_path  <- "inst/extdata/breast/files_metadata.tsv"
downloaded_dir   <- "/vast/scratch/users/shen.m/synapse_data/breast/counts/"
output_dir       <- "/vast/scratch/users/shen.m/htan/hta_2025/0.2.0/counts/"
msk_map_id_path  <- "/vast/scratch/users/shen.m/synapse_data/lung/counts/adata_sample_id_htan_id_map.csv"
# ─────────────────────────────────────────────────────────────────────────────

store <- paste0("~/scratch/htan/", organ, "_all_centers_target_store")

tar_script({
  
  library(targets)
  library(tarchetypes)
  library(tidyverse)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(crew)
  library(crew.cluster)
  
  source("dev/0_center_parser_registry.R")
  
  #organ            <- "breast"          # used only for naming the store
  files_meta_path  <- "/home/users/allstaff/shen.m/projects/HTAN/files_metadata_lung_breastcombined.tsv"
  downloaded_dir   <- "/vast/scratch/users/shen.m/synapse_data/breast_lung_combined/counts/"
  output_dir       <- "/vast/scratch/users/shen.m/htan/hta_2025/0.2.0/counts/"
  msk_map_id_path  <- "/vast/scratch/users/shen.m/synapse_data/breast_lung_combined/counts/adata_sample_id_htan_id_map.csv"
  
  # ── SLURM controllers (small → large with fallback) ─────────────────────
  new_elastic <- function(name, mem_gb, time_min, workers,
                          crashes_max, cpus_per_task = 1, backup = NULL) {
    crew_controller_slurm(
      name             = name,
      workers          = workers,
      crashes_max      = crashes_max,
      seconds_idle     = 30,
      options_cluster  = crew_options_slurm(
        memory_gigabytes_required = mem_gb,
        cpus_per_task             = cpus_per_task,
        time_minutes              = time_min
      ),
      backup = backup
    )
  }
  
  elastic_160 <- new_elastic("elastic_160", 160, 60*24, 8,   2, cpus_per_task = 2)
  elastic_80  <- new_elastic("elastic_80",   80, 60*4,  24,  1, backup = elastic_160)
  elastic_40  <- new_elastic("elastic_40",   40, 60*4,  32,  1, backup = elastic_80)
  elastic_20  <- new_elastic("elastic_20",   20, 60*4,  48,  1, backup = elastic_40)
  elastic_10  <- new_elastic("elastic_10",   10, 60*4,  150, 2, backup = elastic_20)
  elastic_5   <- new_elastic("elastic_5",     5, 60*4,  300, 2, backup = elastic_10)
  
  controllers <- crew_controller_group(
    elastic_5, elastic_10, elastic_20, elastic_40, elastic_80, elastic_160
  )
  
  tar_option_set(
    memory               = "transient",
    garbage_collection   = 100,
    storage              = "worker",
    retrieval            = "worker",
    error                = "continue",
    cue                  = tar_cue(mode = "thorough"),
    format               = "qs",
    workspace_on_error   = TRUE,
    controller           = controllers,
    trust_object_timestamps = TRUE,
    resources = tar_resources(crew = tar_resources_crew(controller = "elastic_5"))
  )
  
  # ── Config (injected by tar_script) ─────────────────────────────────────
  FILES_META_PATH <- files_meta_path
  DOWNLOADED_DIR  <- downloaded_dir
  OUTPUT_DIR      <- output_dir
  MSK_MAP_ID_PATH <- msk_map_id_path
  
  # ── Pipeline ─────────────────────────────────────────────────────────────
  list(
    
    # 1. Read file metadata, attach local paths
    tar_target(
      file_metadata,
      {
        read.csv(FILES_META_PATH, sep = "\t", na.strings = c("NA", ""),
                 header = TRUE) |>
          as_tibble() |>
          mutate(
            Filename_basename = basename(Filename),
            full_path = file.path(DOWNLOADED_DIR, Filename_basename)
          )
      }
    ),
    
    # ══════════════════════════════════════════════════════════════════════
    # HTAN HTAPP  (HTA1) — L3 MTX + L4 TSV
    # ══════════════════════════════════════════════════════════════════════

    tar_target(
      htapp_cell_index_df,
      {
        # L4 TSV batch files → build global cell index
        l4_files <- file_metadata |>
          filter(Atlas.Name == "HTAN HTAPP", Level == "Level 4") |>
          select(sample_id = Biospecimen, processed_file = full_path) |>
          tidyr::separate_rows(sample_id, sep = ",\\s*") |>
          distinct()  

        build_htapp_cell_index(l4_files)
      },
      packages = c("dplyr", "purrr")
    ),

    # Extend htapp_cell_index_df with sequential fallback indices for any
    # sample that has L3 MTX files but no corresponding L4 TSV entry.
    tar_target(
      htapp_cell_index_full_df,
      {
        l4_samples <- unique(htapp_cell_index_df$sample_id)
        missing    <- htapp_metadata |>
          filter(!Biospecimen %in% l4_samples) |>
          select(Biospecimen, barcodes_path) |>
          distinct()
        
        if (nrow(missing) == 0) return(htapp_cell_index_df)
        
        fallback <- purrr::pmap_dfr(missing, function(Biospecimen, barcodes_path) {
          barcodes <- read.delim(barcodes_path, header = FALSE)$V1
          tibble(
            sample_id  = Biospecimen,
            NAME       = barcodes,
            cell_index = seq_along(barcodes)
          )
        })
        
        dplyr::bind_rows(htapp_cell_index_df, fallback)
      },
      packages = c("dplyr", "purrr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),
    
    tar_target(
      htapp_metadata,
      {
        file_metadata |>
          filter(Atlas.Name == "HTAN HTAPP", Level == "Level 3") |>
          mutate(
            file_type = case_when(
              str_detect(Filename_basename, "_matrix\\.mtx\\.gz$")    ~ "mtx_path",
              str_detect(Filename_basename, "_features\\.tsv\\.gz$")  ~ "genes_path",
              str_detect(Filename_basename, "_barcodes\\.tsv\\.gz$")  ~ "barcodes_path"
            ),
            channel_number = Filename_basename |>
              str_replace_all("ch(?=[0-9]+)", "channel") |>
              str_extract("channel[0-9]+")
          ) |>
          group_split(Biospecimen) |>
          map_dfr(~ {
            sg <- .x
            if (nrow(sg) > 4)
              mutate(sg, sample_id = paste(Biospecimen, channel_number, sep = "___"))
            else
              mutate(sg, sample_id = Biospecimen)
          }) |>
          select(full_path, Biospecimen, sample_id, file_type) |>
          pivot_wider(id_cols = c(Biospecimen, sample_id),
                      names_from  = file_type,
                      values_from = full_path)
      },
      packages = c("dplyr", "tidyr", "stringr", "purrr")
    ),

    tar_target(
      htapp_metadata_grouped,
      htapp_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),

    tar_target(
      htapp_sce,
      {
        row      <- htapp_metadata_grouped
        sid      <- row$sample_id[1]
        bid      <- row$Biospecimen[1]
        cid      <- htapp_cell_index_full_df |>
          filter(sample_id == bid) |>
          select(NAME, cell_index)
        parse_HTAPP(sid, row$mtx_path[1], row$genes_path[1],
                    row$barcodes_path[1], cid)
      },
      pattern   = map(htapp_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "Seurat", "dplyr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),

    tar_target(
      save_htapp,
      {
        if (is.null(htapp_sce)) return(NULL)
        save_h5ad(unique(colData(htapp_sce)$sample_id), htapp_sce, OUTPUT_DIR)
      },
      pattern   = map(htapp_sce),
      packages  = c("zellkonverter")
    ),
    
    # ══════════════════════════════════════════════════════════════════════
    # HTAN BU  (HTA3) — L3 CSV + L4 CSV  [Lung only]
    # ══════════════════════════════════════════════════════════════════════
    
    tar_target(
      bu_cell_index_df,
      {
        l4_paths <- file_metadata |>
          filter(Atlas.Name == "HTAN BU") |>
          filter(str_detect(Biospecimen, ",")) |>
          pull(full_path)
        build_bu_cell_index(l4_paths)
      },
      packages = c("dplyr", "purrr")
    ),
    
    tar_target(
      bu_metadata,
      {
        file_metadata |>
          filter(Atlas.Name == "HTAN BU") |>
          filter(str_detect(full_path, "raw")) |>
          select(full_path, sample_id = Biospecimen)
      }
    ),
    
    tar_target(
      bu_metadata_grouped,
      bu_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      bu_sce,
      {
        row <- bu_metadata_grouped
        cid <- bu_cell_index_df |>
          filter(SampleID == row$sample_id[1]) |>
          select(NAME, cell_index)
        parse_BU(row$full_path[1], row$sample_id[1], cid)
      },
      pattern   = map(bu_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "dplyr", "tibble")
    ),
    
    tar_target(
      save_bu,
      {
        if (is.null(bu_sce)) return(NULL)
        save_h5ad(unique(colData(bu_sce)$sample_id), bu_sce, OUTPUT_DIR)
      },
      pattern  = map(bu_sce),
      packages = c("zellkonverter")
    ),
    
    # ══════════════════════════════════════════════════════════════════════
    # HTAN MSK  (HTA8) — L4 h5ad (raw counts embedded)  [Lung only]
    # ══════════════════════════════════════════════════════════════════════
    
    tar_target(
      msk_h5ad_metadata,
      {
        map_id <- read.csv(MSK_MAP_ID_PATH)
        
        h5ad_files <- file_metadata |>
          filter(Atlas.Name == "HTAN MSK", Level == "Level 4",
                 File.Format == "hdf5") |>
          pull(full_path) |> unique()
        
        purrr::map_dfr(h5ad_files, ~ {
          data <- zellkonverter::readH5AD(.x, reader = "R", use_hdf5 = TRUE)
          colData(data) |>
            as.data.frame() |>
            rownames_to_column(".cell") |>
            separate(".cell", c("sample", "barcode"), sep = "_(?=[^_]+$)",
                     remove = FALSE) |>
            left_join(map_id, by = "sample") |>
            dplyr::rename(sample_id = sample_HTAN_ID) |>
            filter(!is.na(sample_id)) |>
            separate_rows(sample_id, sep = ";") |>
            mutate(full_path = .x) |>
            distinct(sample_id, full_path)
        })
      },
      packages = c("dplyr", "purrr", "tidyr", "zellkonverter",
                   "SummarizedExperiment", "tibble", "tidySingleCellExperiment")
    ),
    
    tar_target(
      msk_cell_index_df,
      {
        map_id <- read.csv(MSK_MAP_ID_PATH)
        h5ad_files <- file_metadata |>
          filter(Atlas.Name == "HTAN MSK", Level == "Level 4",
                 File.Format == "hdf5") |>
          pull(full_path) |> unique()
        
        purrr::map_dfr(h5ad_files, ~ {
          data <- zellkonverter::readH5AD(.x, reader = "R", use_hdf5 = TRUE)
          colData(data) |>
            as.data.frame() |>
            rownames_to_column(".cell") |>
            separate(".cell", c("sample", "barcode"), sep = "_(?=[^_]+$)",
                     remove = FALSE) |>
            left_join(map_id, by = "sample") |>
            dplyr::rename(sample_id = sample_HTAN_ID) |>
            filter(!is.na(sample_id)) |>
            separate_rows(sample_id, sep = ";") |>
            mutate(full_path = .x)
        }) |>
          group_by(sample_id) |>
          mutate(cell_index = row_number()) |>
          ungroup() |>
          select(sample_id, full_path, .cell, cell_index)
      },
      packages = c("dplyr", "purrr", "tidyr", "zellkonverter",
                   "SummarizedExperiment", "tibble", "tidySingleCellExperiment")
    ),
    
    tar_target(
      msk_metadata_grouped,
      msk_h5ad_metadata |>
        group_by(sample_id, full_path) |>
        tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      msk_sce_per_file,
      parse_MSK(
        msk_metadata_grouped$full_path,
        msk_metadata_grouped$sample_id,
        msk_cell_index_df |>
          filter(sample_id == msk_metadata_grouped$sample_id,
                 full_path == msk_metadata_grouped$full_path) |>
          select(sample_id, .cell, cell_index),
        MSK_MAP_ID_PATH
      ),
      pattern   = map(msk_metadata_grouped),
      iteration = "list",
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_10")),
      packages  = c("SingleCellExperiment", "zellkonverter", "dplyr", "tidyr", "tibble", "tidySingleCellExperiment")
    ),
    
    tar_target(
      msk_sce_tbl,
      {
        if (is.null(msk_sce_per_file)) return(NULL)
        tibble(sample_id = unique(colData(msk_sce_per_file)$sample_id),
               sce = list(msk_sce_per_file))
      },
      pattern   = map(msk_sce_per_file),
      iteration = "list"
    ),
    
    tar_target(
      msk_sce_by_sample,
      bind_rows(msk_sce_tbl) |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      msk_sce_merged,
      {
        xs           <- msk_sce_by_sample$sce
        common_genes <- Reduce(intersect, map(xs, rownames))
        xs           <- map(xs, ~ {
          s <- .x[common_genes, ]
          SingleCellExperiment(
            assays  = list(counts = assay(s, "counts")),
            colData = colData(s),
            rowData = S4Vectors::DataFrame(row.names = common_genes)
          )
        })
        do.call(cbind, xs)
      },
      pattern   = map(msk_sce_by_sample),
      iteration = "list",
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_10"))
    ),
    
    tar_target(
      save_msk,
      save_h5ad(msk_sce_by_sample$sample_id[1], msk_sce_merged, OUTPUT_DIR),
      pattern  = map(msk_sce_by_sample, msk_sce_merged),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_10")),
      packages = c("zellkonverter")
    )
    ,

    # ══════════════════════════════════════════════════════════════════════
    # HTAN OHSU  (HTA9) — L3 10X HDF5, 1:1 per biospecimen  [Breast]
    # ══════════════════════════════════════════════════════════════════════

    tar_target(
      ohsu_metadata,
      file_metadata |>
        filter(Atlas.Name == "HTAN OHSU", Level == "Level 3",
               File.Format == "hdf5") |>
        select(full_path, sample_id = Biospecimen)
    ),

    tar_target(
      ohsu_metadata_grouped,
      ohsu_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      ohsu_cell_index_df,
      build_ohsu_cell_index(
        ohsu_metadata |> select(sample_id, full_path) |> distinct()
      ),
      packages = c("dplyr", "purrr", "tibble", "rhdf5")
    ),
    
    tar_target(
      ohsu_sce,
      {
        sid <- ohsu_metadata_grouped$sample_id[1]
        cid <- ohsu_cell_index_df |>
          filter(sample_id == sid) |>
          select(NAME, cell_index)
        parse_OHSU(ohsu_metadata_grouped$full_path[1], sid, cid)
      },
      pattern   = map(ohsu_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "DropletUtils", "dplyr")
    ),

    tar_target(
      save_ohsu,
      {
        if (is.null(ohsu_sce)) return(NULL)
        save_h5ad(unique(colData(ohsu_sce)$sample_id), ohsu_sce, OUTPUT_DIR)
      },
      pattern  = map(ohsu_sce),
      packages = c("zellkonverter")
    ),

    # ══════════════════════════════════════════════════════════════════════
    # HTAN WUSTL  (HTA12) — L4 per-sample Seurat RDS  [Breast]
    # L3 MTX not used; RNA assay counts extracted from each RDS directly.
    # Some biospecimens have multiple RDS files — first after filename sort
    # is chosen for reproducibility.
    # ══════════════════════════════════════════════════════════════════════
    
    tar_target(
      wustl_metadata,
      file_metadata |>
        filter(Atlas.Name == "HTAN WUSTL", Level == "Level 4") |>
        select(sample_id = Biospecimen, l4_rds_path = full_path) |>
        group_by(sample_id) |>
        arrange(basename(l4_rds_path), .by_group = TRUE) |>
        slice(1) |>
        ungroup(),
      packages = c("dplyr")
    ),
    
    tar_target(
      wustl_metadata_grouped,
      wustl_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      wustl_cell_index_df,
      build_wustl_cell_index(
        wustl_metadata |> select(sample_id, l4_rds_path)
      ),
      pattern = map(wustl_metadata_grouped),
      packages = c("Seurat", "dplyr", "purrr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_40"))
    ),
    
    tar_target(
      wustl_sce,
      {
        sid <- wustl_metadata_grouped$sample_id[1]
        cid <- wustl_cell_index_df |>
          filter(sample_id == sid) |>
          select(NAME, cell_index)
        parse_WUSTL(sid, wustl_metadata_grouped$l4_rds_path[1], cid)
      },
      pattern   = map(wustl_metadata_grouped, wustl_cell_index_df),
      iteration = "list",
      packages  = c("SingleCellExperiment", "Seurat", "dplyr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_40"))
    ),
    
    tar_target(
      save_wustl,
      {
        if (is.null(wustl_sce)) return(NULL)
        save_h5ad(unique(colData(wustl_sce)$sample_id), wustl_sce, OUTPUT_DIR)
      },
      pattern  = map(wustl_sce),
      packages = c("zellkonverter")
    ),
    # 
    # ══════════════════════════════════════════════════════════════════════
    # HTAN Duke  (HTA6) — L4 merged RDS (22 samples) + L3 MTX fallback
    #                     for the 8 biospecimens absent from L4
    # ══════════════════════════════════════════════════════════════════════
    
    tar_target(
      duke_l4_rds_path,
      file_metadata |>
        filter(Atlas.Name == "HTAN Duke", Level == "Level 4",
               File.Format == "RData") |>
        pull(full_path) |> unique() |> (\(x) x[1])()
    ),
    
    tar_target(
      duke_metadata,
      {
        file_metadata |>
          filter(Atlas.Name == "HTAN Duke", Level == "Level 3") |>
          mutate(file_type = case_when(
            str_detect(Filename_basename, "-matrix\\.mtx\\.gz$")    ~ "mtx_path",
            str_detect(Filename_basename, "-features\\.tsv\\.gz$")  ~ "genes_path",
            str_detect(Filename_basename, "-barcodes\\.tsv\\.gz$")  ~ "barcodes_path"
          )) |>
          select(full_path, sample_id = Biospecimen, file_type) |>
          pivot_wider(id_cols = sample_id,
                      names_from  = file_type,
                      values_from = full_path)
      },
      packages = c("dplyr", "tidyr", "stringr")
    ),
    
    # Unified cell index: L4 Seurat entries (sorted barcodes → reproducible
    # cell_index) + sequential fallback entries for L3-only samples.
    tar_target(
      duke_cell_index_df,
      build_duke_cell_index(duke_l4_rds_path, duke_metadata),
      packages = c("Seurat", "dplyr", "purrr", "tibble")
    ),
    
    # ── L4 path: parse all 22 Seurat-covered samples in one shot ─────────
    # Returns a named list of SCEs (name = sample_id). Loading the RDS once
    # avoids redundant IO vs. per-sample branching.
    tar_target(
      duke_l4_sce,
      parse_Duke_seurat(duke_l4_rds_path, duke_cell_index_df),
      packages  = c("Seurat", "SingleCellExperiment", "dplyr", "purrr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_80"))
    ),
    
    tar_target(
      save_duke_l4,
      {
        if (is.null(duke_l4_sce)) return(NULL)
        save_h5ad(unique(colData(duke_l4_sce[[1]])$sample_id), duke_l4_sce[[1]], OUTPUT_DIR)
      },
      pattern = map(duke_l4_sce),
      packages = c("zellkonverter", "purrr")
    ),
    
    # ── L3-only fallback: samples in L3 metadata absent from the L4 RDS ──
    tar_target(
      duke_l3_only_metadata,
      {
        l4_samples <- unique(duke_cell_index_df$sample_id[
          duke_cell_index_df$sample_id %in%
            names(duke_l4_sce)[!sapply(duke_l4_sce, is.null)]
        ])
        duke_metadata |> filter(!sample_id %in% l4_samples)
      },
      packages = c("dplyr")
    ),
    
    tar_target(
      duke_l3_only_metadata_grouped,
      duke_l3_only_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      duke_l3_sce,
      {
        row <- duke_l3_only_metadata_grouped
        sid <- row$sample_id[1]
        cid <- duke_cell_index_df |>
          filter(sample_id == sid) |>
          select(NAME, cell_index)
        parse_Duke(sid, row$mtx_path[1], row$genes_path[1],
                   row$barcodes_path[1], cid)
      },
      pattern   = map(duke_l3_only_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "Seurat", "dplyr", "tibble")
    ),
    
    tar_target(
      save_duke_l3,
      {
        if (is.null(duke_l3_sce)) return(NULL)
        save_h5ad(unique(colData(duke_l3_sce)$sample_id), duke_l3_sce, OUTPUT_DIR)
      },
      pattern  = map(duke_l3_sce),
      packages = c("zellkonverter")
    ),
    # 
    # ══════════════════════════════════════════════════════════════════════
    # HTAN DFCI  (HTA5) — POOLED, no cell-level demultiplexing
    # ══════════════════════════════════════════════════════════════════════

    tar_target(
      dfci_metadata,
      file_metadata |>
        filter(Atlas.Name == "HTAN DFCI", Level == "Level 3") |>
        select(full_path, #pool_id = Filename_basename, 
               sample_id = Biospecimen) |> 
        separate_rows(sample_id, sep = ",\\s*")
    ),
    tar_target(
      dfci_metadata_grouped,
      dfci_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),

    tar_target(
      dfci_cell_index_df,
      build_dfci_cell_index(
        dfci_metadata |> select(sample_id, full_path) |> distinct()
      ),
      packages = c("dplyr", "purrr", "tibble", "rhdf5")
    ),

    tar_target(
      dfci_sce,
      {
        sid <- dfci_metadata_grouped$sample_id[1]
        cid <- dfci_cell_index_df |>
          filter(sample_id == sid) |>
          select(NAME, cell_index)
        parse_DFCI(dfci_metadata_grouped$full_path[1], sid, cid)
      },
      pattern   = map(dfci_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "DropletUtils", "dplyr")
    ),

    tar_target(
      save_dfci,
      {
        if (is.null(dfci_sce)) return(NULL)
        save_h5ad(unique(colData(dfci_sce)$sample_id), dfci_sce, OUTPUT_DIR)
      },
      pattern  = map(dfci_sce),
      packages = c("zellkonverter")
    )

  )
  
}, script = paste0(store, "_target_script.R"), ask = FALSE)


job::job({
  tar_make(
    script   = paste0(store, "_target_script.R"),
    store    = store,
    reporter = "summary"
  )
})

tar_workspace(save_duke_l4_bc70e216b71c6929,store = store, script   = paste0(store, "_target_script.R"))
tar_progress_branches(store = store) |> mutate(pending = branches - skipped - completed)
