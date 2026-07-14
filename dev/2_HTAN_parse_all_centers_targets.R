# =============================================================================
# Unified HTAN targets pipeline — dispatches per center using the registry
#
# Usage:
#   Set the metadata paths below, then source this file.
#   Downloaded files are organised by center_id subdirectory under
#   `downloaded_dir` (matching the layout produced by script 1).
#
# Outputs: one .h5ad per HTAN biospecimen ID, saved to `output_dir`.
# =============================================================================

library(targets)
library(dplyr)

# # Bind files_metadata
# read.csv("/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv",
#          sep = "\t", na.strings = c("NA",""), header = TRUE) |> bind_rows(
#            read.csv("/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/breast/files_metadata.tsv",
#                     sep = "\t", na.strings = c("NA",""), header = TRUE) 
#          ) |> as_tibble() |>
#   write.table(file = "/home/users/allstaff/shen.m/projects/HTAN/files_metadata_lung_breastcombined.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = "")


# ── Configure these lines per run ────────────────────────────────────────────
organ            <- "all_centers"     # used only for naming the store
files_meta_path  <- "inst/extdata/files_metadata_scRNAseq_synapse_level3_4.tsv"
downloaded_dir   <- "/vast/scratch/users/shen.m/synapse_data/"   # center_id subdirs live here
output_dir       <- "/vast/scratch/users/shen.m/htan/hta_2025/0.3.0/parsed_counts/"
msk_map_id_path  <- "/vast/scratch/users/shen.m/synapse_data/breast_lung_combined/counts/adata_sample_id_htan_id_map.csv"
# ─────────────────────────────────────────────────────────────────────────────

store <- paste0("/vast/scratch/users/shen.m/htan/", organ, "_all_centers_target_store")

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
  files_meta_path  <- "/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/files_metadata_scRNAseq_synapse_level3_4.tsv"
  downloaded_dir   <- "/vast/scratch/users/shen.m/synapse_data/"
  output_dir       <- "/vast/scratch/users/shen.m/htan/hta_2025/0.3.0/parsed_counts/"
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
    
    # 1. Read file metadata, attach local paths.
    #    Most files: DOWNLOADED_DIR/<center_id>/<basename>  (flat per center).
    #    CHOP L3 MTX exception: DOWNLOADED_DIR/HTAN_CHOP/<run_dir>___<basename>
    #      — script 1 renames each file as  <run_dir>___<basename>  (flat, no
    #        subdirectory) so that barcodes/features/matrix from different run
    #        directories never collide in HTAN_CHOP/.
    #        Example: Other_2263_CD33_Neg_scRNA___barcodes.tsv.gz
    tar_target(
      file_metadata,
      {
        chop_l3_mtx <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
        read.csv(FILES_META_PATH, sep = "\t", na.strings = c("NA", ""),
                 header = TRUE) |>
          as_tibble() |>
          mutate(
            Filename_basename = basename(Filename),
            center_id         = sub(" ", "_", Atlas.Name),
            full_path         = if_else(
              Atlas.Name == "HTAN CHOP" & Level == "Level 3" &
                Filename_basename %in% chop_l3_mtx,
              file.path(DOWNLOADED_DIR, center_id,
                        paste(basename(dirname(Filename)), Filename_basename,
                              sep = "___")),
              file.path(DOWNLOADED_DIR, center_id, Filename_basename)
            )
          )
      }
    ),
    

    # ══════════════════════════════════════════════════════════════════════
    # HTAN HTAPP  (HTA1) — L3 MTX + L4 TSV
    # ══════════════════════════════════════════════════════════════════════
    tar_target(
      htapp_metadata,
      {
        file_metadata |>
          filter(Atlas.Name == "HTAN HTAPP", Level == "Level 3") |>
          mutate(
            file_type = case_when(
              str_detect(Filename_basename, "_matrix\\.mtx(\\.gz)?$")              ~ "mtx_path",
              str_detect(Filename_basename, "_(features|genes)\\.tsv(\\.gz)?$")    ~ "genes_path",
              str_detect(Filename_basename, "_barcodes\\.tsv(\\.gz)?$")            ~ "barcodes_path"
            ),
            filename_renamed = Filename_basename |>
              str_replace_all("ch(?=[0-9]+)", "channel"),
            channel_number = filename_renamed |>
              str_extract("channel[0-9]+"),
            library_id = filename_renamed |>
              str_extract("[^_]+(?=_channel[0-9]+)"),
            library_channel = str_c(library_id, "_", channel_number)
          ) |>
          group_split(Biospecimen) |>
          map_dfr(~ {
            sg <- .x
            if (nrow(sg) > 4)
              mutate(sg, sample_id = paste(Biospecimen, library_channel, sep = "___"))
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
          select(NAME, cell_index, sample_id)
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
          filter(SampleID == row$sample_id[1])
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
        filter(Atlas.Name == "HTAN OHSU", Level == "Level 3", # Use Level 3 here because Level 4 anndata doesnt hold meaningful colData
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
          filter(sample_id == sid)
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
          filter(sample_id == sid)
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
          filter(sample_id == sid)
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
          filter(sample_id == sid)
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
    ),
    # ══════════════════════════════════════════════════════════════════════
    # HTAN CHOP  (HTA4) — L4 per-sample RDS (NBL) + L3 MTX fallback (ped-glioma)
    #
    # L4: seurat_objects_count_based/seurat_object_NB_7767_*.rds
    #     Single biospecimen per file; integrated/merged RDS files are ignored.
    # L3: Standard CellRanger MTX triplet (barcodes/features/matrix) — used for
    #     all samples not covered by a per-sample L4 RDS.
    # ══════════════════════════════════════════════════════════════════════
    
    tar_target(
      chop_l4_rds_metadata,
      file_metadata |>
        filter(Atlas.Name == "HTAN CHOP", Level == "Level 4",
               File.Format == "RData",
               !grepl(",", Biospecimen),
               Biospecimen != "0 Biospecimens",
               str_detect(Filename, "seurat_objects_count_based")) |>
        select(sample_id = Biospecimen, l4_rds_path = full_path),
      packages = c("dplyr", "stringr")
    ),
    
    tar_target(
      chop_l3_metadata,
      {
        # snRNA directories (ped-glioma snRNA-seq) are excluded entirely.
        # Each biospecimen+run_dir combination that forms a complete triplet
        # is treated as an independent sample.
        # sample_id: plain Biospecimen for single-run samples;
        #            paste(Biospecimen, run_dir, sep="___") for multi-run
        #            (e.g. HSPC donors with CD34/LinCD38CD34/Live sorts).
        file_metadata |>
          filter(Atlas.Name == "HTAN CHOP", Level == "Level 3",
                 !grepl(",", Biospecimen),
                 Biospecimen != "0 Biospecimens",
                 !str_detect(Filename, "_bak/"),
                 !str_detect(Filename, "snRNA")) |>
          mutate(
            file_type = case_when(
              Filename_basename == "barcodes.tsv.gz" ~ "barcodes_path",
              Filename_basename == "features.tsv.gz" ~ "genes_path",
              Filename_basename == "matrix.mtx.gz"   ~ "mtx_path"
            ),
            run_dir = basename(dirname(Filename))
          ) |>
          filter(!is.na(file_type)) |>
          select(full_path, biospecimen = Biospecimen, run_dir, file_type) |>
          group_by(biospecimen, run_dir, file_type) |>
          filter(n() == 1L) |>       # drop any incomplete triplet
          ungroup() |>
          pivot_wider(id_cols     = c(biospecimen, run_dir),
                      names_from  = file_type,
                      values_from = full_path) |>
          filter(!is.na(barcodes_path) & !is.na(genes_path) & !is.na(mtx_path)) |>
          group_by(biospecimen) |>
          mutate(sample_id = if (n() > 1L)
            paste(biospecimen, run_dir, sep = "___")
            else
              biospecimen) |>
          ungroup() |>
          select(sample_id, biospecimen, run_dir, barcodes_path, genes_path, mtx_path)
      },
      packages = c("dplyr", "tidyr", "stringr")
    ),
    
    tar_target(
      chop_l4_rds_metadata_grouped,
      chop_l4_rds_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    # ── Cell index: L4 per-sample (one RDS per branch) ───────────────────
    tar_target(
      chop_l4_cell_index_df,
      build_chop_l4_cell_index(
        chop_l4_rds_metadata_grouped |> select(sample_id, l4_rds_path)
      ),
      pattern   = map(chop_l4_rds_metadata_grouped),
      packages  = c("Seurat", "dplyr", "purrr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),
    
    # ── L4 path: parse per-sample RDS ────────────────────────────────────
    tar_target(
      chop_l4_sce,
      {
        sid <- chop_l4_rds_metadata_grouped$sample_id[1]
        parse_CHOP_seurat(sid, chop_l4_rds_metadata_grouped$l4_rds_path[1],
                          chop_l4_cell_index_df)
      },
      pattern   = map(chop_l4_rds_metadata_grouped, chop_l4_cell_index_df),
      iteration = "list",
      packages  = c("SingleCellExperiment", "Seurat", "dplyr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),
    
    tar_target(
      save_chop_l4,
      {
        if (is.null(chop_l4_sce)) return(NULL)
        save_h5ad(unique(colData(chop_l4_sce)$sample_id), chop_l4_sce, OUTPUT_DIR)
      },
      pattern  = map(chop_l4_sce),
      packages = c("zellkonverter")
    ),
    
    # ── L3 fallback: samples absent from per-sample L4 RDS ───────────────
    # Compare on biospecimen (plain HTAN ID) because L4 metadata uses plain
    # IDs while L3 sample_id may carry a run_dir suffix for multi-run samples.
    tar_target(
      chop_l3_only_metadata,
      chop_l3_metadata |>
        filter(!biospecimen %in% chop_l4_rds_metadata$sample_id),
      packages = c("dplyr")
    ),
    
    tar_target(
      chop_l3_cell_index_df,
      build_chop_l3_cell_index(
        chop_l3_only_metadata |> select(sample_id, barcodes_path) |> distinct()
      ),
      packages = c("dplyr", "purrr", "tibble")
    ),
    
    tar_target(
      chop_l3_only_metadata_grouped,
      chop_l3_only_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      chop_l3_sce,
      {
        row <- chop_l3_only_metadata_grouped
        sid <- row$sample_id[1]
        cid <- chop_l3_cell_index_df |>
          filter(sample_id == sid) |>
          select(NAME, cell_index, sample_id)
        parse_CHOP(sid, row$mtx_path[1], row$genes_path[1],
                   row$barcodes_path[1], cid)
      },
      pattern   = map(chop_l3_only_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "Seurat", "dplyr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),
    
    tar_target(
      save_chop_l3,
      {
        if (is.null(chop_l3_sce)) return(NULL)
        save_h5ad(unique(colData(chop_l3_sce)$sample_id), chop_l3_sce, OUTPUT_DIR)
      },
      pattern  = map(chop_l3_sce),
      packages = c("zellkonverter")
    ),
    # ══════════════════════════════════════════════════════════════════════
    # HTAN Stanford  (HTA10) — L3 MTX triplet, no L4  [Colon NOS]
    # "-R0" re-run files share the same Biospecimen ID as the primary run
    # and are excluded so each biospecimen maps to exactly one MTX triplet.
    # ══════════════════════════════════════════════════════════════════════
    
    tar_target(
      stanford_metadata,
      {
        file_metadata |>
          filter(Atlas.Name == "HTAN Stanford", Level == "Level 3",
                 # R0 is a replicate run
                 !str_detect(Filename_basename, "R0")) |>
          mutate(file_type = case_when(
            str_detect(Filename_basename, "_matrix\\.mtx\\.gz$")    ~ "mtx_path",
            str_detect(Filename_basename, "_features\\.tsv\\.gz$")  ~ "genes_path",
            str_detect(Filename_basename, "_barcodes\\.tsv\\.gz$")  ~ "barcodes_path"
          )) |>
          filter(!is.na(file_type)) |>
          select(full_path, sample_id = Biospecimen, file_type) |>
          pivot_wider(id_cols     = sample_id,
                      names_from  = file_type,
                      values_from = full_path) |>
          filter(!is.na(barcodes_path) & !is.na(genes_path) & !is.na(mtx_path))
      },
      packages = c("dplyr", "tidyr", "stringr")
    ),
    
    tar_target(
      stanford_metadata_grouped,
      stanford_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      stanford_cell_index_df,
      build_stanford_cell_index(
        stanford_metadata |> select(sample_id, barcodes_path) |> distinct()
      ),
      packages = c("dplyr", "purrr", "tibble")
    ),
    
    tar_target(
      stanford_sce,
      {
        row <- stanford_metadata_grouped
        sid <- row$sample_id[1]
        cid <- stanford_cell_index_df |>
          filter(sample_id == sid) |>
          select(NAME, cell_index, sample_id)
        parse_stanford(sid, row$mtx_path[1], row$genes_path[1],
                       row$barcodes_path[1], cid)
      },
      pattern   = map(stanford_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "Seurat", "DropletUtils", "dplyr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),
    
    tar_target(
      save_stanford,
      {
        if (is.null(stanford_sce)) return(NULL)
        save_h5ad(unique(colData(stanford_sce)$sample_id), stanford_sce, OUTPUT_DIR)
      },
      pattern  = map(stanford_sce),
      packages = c("zellkonverter")
    ),
    
    # ══════════════════════════════════════════════════════════════════════
    # HTAN Vanderbilt  (HTA11) — L3 CSV, cell-type partitioned  [Colon]
    #
    # CSV path A: YX series — always 2 files per biospecimen (_epi + _immune),
    #             cbind-ed into one SCE.
    # CSV path B: as series — 1–4 CSVs from different study batches sharing
    #             the same Biospecimen ID, cbind-ed after gene-set intersection.
    # MTX path A: P9142 (2 biospecimens, single IDs) — parsed via parse_stanford.
    # MTX path B: Lau series (multi-ID pooled lanes) — POOLED, flagged/skipped.
    # ══════════════════════════════════════════════════════════════════════
    
    tar_target(
      vanderbilt_csv_metadata,
      {
        file_metadata |>
          filter(Atlas.Name == "HTAN Vanderbilt", Level == "Level 3",
                 File.Format == "csv") |>
          select(sample_id = Biospecimen, csv_path = full_path) |>
          group_by(sample_id) |>
          summarise(csv_paths = list(sort(csv_path)), .groups = "drop")
      },
      packages = c("dplyr")
    ),
    
    tar_target(
      vanderbilt_csv_metadata_grouped,
      vanderbilt_csv_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      vanderbilt_csv_cell_index_df,
      build_vanderbilt_csv_cell_index(
        vanderbilt_csv_metadata |> select(sample_id, csv_paths)
      ),
      packages = c("dplyr", "purrr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),
    
    tar_target(
      vanderbilt_csv_sce,
      {
        row <- vanderbilt_csv_metadata_grouped
        sid <- row$sample_id[1]
        cid <- vanderbilt_csv_cell_index_df |>
          filter(sample_id == sid) |>
          select(NAME, cell_index, sample_id)
        parse_vanderbilt_csv(sid, row$csv_paths[[1]], cid)
      },
      pattern   = map(vanderbilt_csv_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "dplyr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),
    
    tar_target(
      save_vanderbilt_csv,
      {
        if (is.null(vanderbilt_csv_sce)) return(NULL)
        save_h5ad(unique(colData(vanderbilt_csv_sce)$sample_id),
                  vanderbilt_csv_sce, OUTPUT_DIR)
      },
      pattern  = map(vanderbilt_csv_sce),
      packages = c("zellkonverter")
    ),
    
    # ── P9142 MTX: 2 single-biospecimen samples, parsed like Stanford ────
    tar_target(
      vanderbilt_mtx_metadata,
      {
        file_metadata |>
          filter(Atlas.Name == "HTAN Vanderbilt", Level == "Level 3",
                 !grepl(",", Biospecimen),               # exclude pooled Lau lanes
                 File.Format %in% c("mtx", "tsv")) |>
          mutate(file_type = case_when(
            str_detect(Filename_basename, "_matrix\\.mtx\\.gz$")    ~ "mtx_path",
            str_detect(Filename_basename, "_features\\.tsv\\.gz$")  ~ "genes_path",
            str_detect(Filename_basename, "_barcodes\\.tsv\\.gz$")  ~ "barcodes_path"
          )) |>
          filter(!is.na(file_type)) |>
          select(full_path, sample_id = Biospecimen, file_type) |>
          pivot_wider(id_cols     = sample_id,
                      names_from  = file_type,
                      values_from = full_path) |>
          filter(!is.na(barcodes_path) & !is.na(genes_path) & !is.na(mtx_path))
      },
      packages = c("dplyr", "tidyr", "stringr")
    ),
    
    tar_target(
      vanderbilt_mtx_cell_index_df,
      build_stanford_cell_index(
        vanderbilt_mtx_metadata |> select(sample_id, barcodes_path)
      ),
      packages = c("dplyr", "purrr", "tibble")
    ),
    
    tar_target(
      vanderbilt_mtx_metadata_grouped,
      vanderbilt_mtx_metadata |> group_by(sample_id) |> tar_group(),
      iteration = "group"
    ),
    
    tar_target(
      vanderbilt_mtx_sce,
      {
        row <- vanderbilt_mtx_metadata_grouped
        sid <- row$sample_id[1]
        cid <- vanderbilt_mtx_cell_index_df |>
          filter(sample_id == sid) |>
          select(NAME, cell_index, sample_id)
        parse_vanderbilt(sid, row$mtx_path[1], row$genes_path[1],
                         row$barcodes_path[1], cid)
      },
      pattern   = map(vanderbilt_mtx_metadata_grouped),
      iteration = "list",
      packages  = c("SingleCellExperiment", "Seurat", "DropletUtils", "dplyr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ),
    
    tar_target(
      save_vanderbilt_mtx,
      {
        if (is.null(vanderbilt_mtx_sce)) return(NULL)
        save_h5ad(unique(colData(vanderbilt_mtx_sce)$sample_id),
                  vanderbilt_mtx_sce, OUTPUT_DIR)
      },
      pattern  = map(vanderbilt_mtx_sce),
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

tar_workspace(vanderbilt_csv_sce_8843e4ca6708c8b3, store = store, script   = paste0(store, "_target_script.R"))
tar_progress_branches(store = store) |> mutate(pending = branches - skipped - completed) |>
  print(n=50)


