# =============================================================================
# HTAN Multi-Organ Harmonisation — Center-Based Parser Registry
# =============================================================================
#
# CORE INSIGHT
# ------------
# HTAN data heterogeneity is BY CENTER, not by organ. The same center submits
# the same file format for Lung, Breast, Colon, etc. Adding a new organ
# therefore only requires identifying which centers contributed and checking
# if their format has already been handled.
#
# DECISION TREE
# -------------
# For each center, ask two questions:
#   1. Does the Level 4 file contain raw counts?
#      YES → read from L4 (skip L3 entirely)         [MSK: L4 h5ad]
#      NO  → read counts from L3, use L4 as a filter  [HTAPP, BU, Duke, WUSTL]
#      L4 absent → read all cells from L3             [OHSU, DFCI*]
#
#   2. If reading from L3, what format?
#      MTX triplet  → ReadMtx / read10xCounts         [HTAPP, Duke, WUSTL]
#      Dense CSV    → read.csv                         [BU, some MSK]
#      10X HDF5     → DropletUtils::read10xCounts      [OHSU, DFCI]
#
# * DFCI is a POOLED 10X experiment: multiple HTAN biospecimen IDs share one
#   HDF5 file with no per-cell sample assignment available in the metadata.
#   Cannot be demultiplexed here; flagged as pool — needs separate handling.
#
# CENTER SUMMARY (Lung + Breast + Paediatric as of 2025-10)
# ----------------------------------------------------------
#  Center      | Atlas ID | Organ(s)             | L3 format        | L4 format           | Strategy
#  ------------|----------|----------------------|------------------|---------------------|-----------------------------
#  HTAN HTAPP  | HTA1     | Lung, Breast         | MTX triplet      | TSV (cell list)     | L3_MTX + L4_TSV
#  HTAN BU     | HTA3     | Lung                 | Dense CSV        | CSV (cell list)     | L3_CSV + L4_CSV
#  HTAN CHOP   | HTA4     | NBL, ped-glioma      | MTX triplet      | RDS (per smp, NBL)  | L4_RDS + L3_MTX_fallback
#  HTAN MSK    | HTA8     | Lung                 | MTX / CSV        | h5ad (counts)       | L4_H5AD  ← preferred
#  HTAN DFCI   | HTA5     | Breast               | 10X HDF5 (pool)  | none                | POOLED — flag
#  HTAN Duke   | HTA6     | Breast               | MTX triplet      | RDS (all smpls)     | L3_MTX + L4_RDS
#  HTAN OHSU   | HTA9     | Breast               | 10X HDF5 (1:1)   | none                | L3_10X_HDF5
#  HTAN WUSTL  | HTA12    | Breast               | MTX triplet      | RDS (per sample)    | L3_MTX + L4_RDS
#
# WHY L4_H5AD FOR MSK?
# ---------------------
# MSK L4 h5ads are merged, QC-filtered AnnData objects produced by the
# submitting lab (Laughney/Bhatt, Chan, Glasner cohorts). They hold raw counts
# in assay[1] and cover 121 biospecimens — including 3 that have NO L3 file.
# Reading one h5ad per cohort is far cheaper than reading 110-118 MTX/CSV
# files and then cross-referencing cell IDs. The HTAPP "Unknown" h5ads in
# Breast are the OPPOSITE: batch-corrected integration outputs (not raw counts)
# and a strict subset of the MTX triplets — ignore them.
#
# HOW TO ADD A NEW ORGAN
# ----------------------
# 1. Download files_metadata.tsv for the new organ from data.humantumoratlas.org
# 2. Run: table(files_metadata$Atlas.Name, files_metadata$Level, files_metadata$File.Format)
# 3. Compare Atlas.Name values against CENTRE_REGISTRY below
#    - Already registered → zero new code needed
#    - New center        → write parse_<CENTER>() following the interface below
#                          and add an entry to CENTRE_REGISTRY
# 4. Run the unified targets pipeline pointing at the new metadata file
#
# STANDARD PARSER INTERFACE
# --------------------------
# Every parse_* function returns a SingleCellExperiment with:
#   - assayNames: "counts"  (raw integer counts)
#   - rownames:   Ensembl gene IDs (ENSG...)
#   - colnames:   numeric cell index (integer, assigned globally per sample_id)
#   - colData:    at minimum $sample_id (HTAN biospecimen ID)
# =============================================================================

library(SingleCellExperiment)
library(zellkonverter)
library(BiocParallel)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)


# -----------------------------------------------------------------------------
# Shared utilities
# -----------------------------------------------------------------------------

#' Convert gene symbols in rownames to Ensembl IDs using EnsDb.Hsapiens.v86
convert_gene_to_ensembl <- function(sce) {
  library(AnnotationFilter)
  gene_id <- ensembldb::mapIds(
    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    keys      = rownames(sce),
    column    = "GENEID",
    keytype   = "GENENAME",
    multiVals = "first"
  )
  keep <- !is.na(gene_id)
  sce  <- sce[keep, ]
  rowData(sce)$gene_id <- gene_id[keep]
  rownames(sce)        <- gene_id[keep]
  sce
}

#' Assign numeric cell indices within a sample, rename colnames, and ensure
#' rownames are Ensembl IDs. Expects colData to already contain $sample_id.
#'
#' @param sce        SingleCellExperiment
#' @param cell_ids   Named character vector: names = old barcode, values = new index.
#'                   If NULL, indices are assigned by row_number().
finalise_sce <- function(sce, cell_ids = NULL) {
  if (!is.null(cell_ids)) {
    sce <- sce[, names(cell_ids)]
    colnames(sce) <- as.character(cell_ids)
  }
  if (!all(grepl("^ENSG", rownames(sce)))) {
    sce <- convert_gene_to_ensembl(sce)
  }
  # strip reduced dims and metadata to keep the object lean
  reducedDims(sce) <- list()
  metadata(sce)    <- list()
  sce
}

# OOM-signature patterns from killed BiocParallel children.
# Extend this list as new variants appear in HPC logs — better to escalate
# too eagerly on small tiers than to silently error on a genuine OOM.
.EMPTYDROP_OOM_PATTERNS <- c(
  "wrong args for environment subassignment"
)

#' Filter empty droplets from a raw SCE (MTX triplet or 10X HDF5 input).
#' Applied before any cell-index join so that the cell set passed downstream
#' contains only real cells.  FDR NAs (very low-count barcodes) are treated as
#' non-cells, matching DropletUtils documentation.
#'
#' Error handling:
#'   - "no counts available …" → all barcodes are real cells; skip filter.
#'   - OOM signatures         → exit with status 137 so the targets/mirai
#'                               backend registers a crash and escalates to a
#'                               higher-memory HPC tier.
#'   - Any other error        → re-thrown as a normal target failure.
filter_empty_droplets <- function(sce, workers = 5) {
  counts(sce) <- as(counts(sce), "CsparseMatrix")
  set.seed(100)
  e.out <- tryCatch(
    DropletUtils::emptyDrops(
      counts(sce),
      BPPARAM = MulticoreParam(workers = workers)
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      
      if (grepl("no counts available to estimate the ambient profile", msg,
                fixed = TRUE)) {
        message("emptyDrops: all barcodes exceed the ambient threshold — skipping filter.")
        return(NULL)
      }
      
      is_oom <- any(vapply(.EMPTYDROP_OOM_PATTERNS, grepl, logical(1),
                           x = msg, ignore.case = TRUE))
      if (is_oom) {
        message("emptyDrops: likely OOM in BiocParallel child — exiting with ",
                "status 137 to trigger mirai resource escalation.")
        quit(save = "no", status = 137L, runLast = FALSE)
      }
      
      stop(e)
    }
  )
  if (is.null(e.out)) return(sce)
  is.cell <- !is.na(e.out$FDR) & e.out$FDR <= 0.01
  sce[, which(is.cell)]
}

#' Strip numeric version suffixes (.1, .2, …) from a gene-name vector.
#' Suffixes are removed ONLY when the result remains unique across all genes.
#' If stripping would introduce duplicates (i.e. the suffixes are meaningful
#' rather than CellRanger disambiguation artefacts), the original vector is
#' returned unchanged with a warning.
strip_gene_version_suffix <- function(genes) {
  stripped <- sub("\\.[0-9]+$", "", genes)
  if (!anyDuplicated(stripped)) return(stripped)
  warning("strip_gene_version_suffix: stripping would produce ",
          sum(duplicated(stripped)), " duplicate gene names — keeping originals.")
  genes
}

#' Save a single SCE as a gzip-compressed h5ad
save_h5ad <- function(sample_id, sce, save_directory) {
  if (!dir.exists(save_directory)) dir.create(save_directory, recursive = TRUE)
  out <- file.path(save_directory, paste0(sample_id, ".h5ad"))
  zellkonverter::writeH5AD(sce, file = out, compression = "gzip")
  message("saved: ", out)
  out
}

#' Resolve a file path case-insensitively on case-sensitive filesystems (HPC).
#' Synapse sometimes downloads files with a different case than the metadata
#' records (e.g. metadata says "4563as25.csv", disk has "4563AS25.csv").
#' Returns the path unchanged if it already exists; otherwise scans the parent
#' directory for a case-insensitive basename match and returns the first hit.
#' Stops with an informative message if no match is found at all.
resolve_csv_path <- function(path) {
  if (file.exists(path)) return(path)
  dir  <- dirname(path)
  base <- basename(path)
  # Case-insensitive match by comparing lowercased basenames directly —
  # avoids regex metacharacter escaping issues entirely.
  all_files <- list.files(dir, full.names = TRUE)
  hits <- all_files[tolower(basename(all_files)) == tolower(base)]
  if (length(hits) == 1L) {
    message("resolve_csv_path: '", base, "' not found; using '",
            basename(hits), "' (case-insensitive match).")
    return(hits)
  }
  if (length(hits) > 1L) {
    warning("resolve_csv_path: multiple case-insensitive matches for '", base,
            "' in ", dir, "; using first: ", basename(hits[1]))
    return(hits[1])
  }
  stop("resolve_csv_path: file not found and no case-insensitive match: ", path)
}

# -----------------------------------------------------------------------------
# Parser: HTAN HTAPP  (HTA1)
# Strategy: L3 MTX triplet  +  L4 TSV cell filter
#
# L3: CellRanger MTX triplet per channel.
#     Multi-channel biospecimens → sample_id = paste(Biospecimen, channel, sep="___")
# L4: TSV with columns NAME (barcode), Biospecimen, HTAPP_id, …
#     One L4 TSV covers multiple biospecimens (batch file).
#     The NAME column encodes the barcode as output by CellRanger.
# -----------------------------------------------------------------------------

# Reads the mtx header to get expected dimensions, without loading the whole matrix
get_mtx_dims <- function(mtx_path) {
  con <- if (grepl("\\.gz$", mtx_path)) gzfile(mtx_path, "r") else file(mtx_path, "r")
  on.exit(close(con))
  line <- readLines(con, n = 1)
  while (startsWith(line, "%")) line <- readLines(con, n = 1)
  dims <- as.integer(strsplit(trimws(line), "\\s+")[[1]])
  list(n_features = dims[1], n_barcodes = dims[2])
}

# Wraps ReadMtx: tries as-is, and on a row/col count mismatch, measures the
# actual line-count offset per file and retries with the correct skip values
read_mtx_safe <- function(mtx_path, genes_path, barcodes_path) {
  tryCatch({
    Seurat::ReadMtx(mtx = mtx_path, cells = barcodes_path,
                    features = genes_path, feature.column = 1)
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (!grepl("but found", msg)) stop(e)  # not the mismatch we know how to fix
    
    dims <- get_mtx_dims(mtx_path)
    n_genes    <- length(readLines(genes_path))
    n_barcodes <- length(readLines(barcodes_path))
    
    skip.feature <- n_genes    - dims$n_features
    skip.cell    <- n_barcodes - dims$n_barcodes
    
    if (skip.feature < 0 || skip.cell < 0) {
      stop(sprintf(
        "Unresolvable mismatch for %s: mtx expects %d features/%d barcodes, found %d/%d lines (negative skip)",
        basename(mtx_path), dims$n_features, dims$n_barcodes, n_genes, n_barcodes
      ))
    }
    
    message(sprintf(
      "%s: retrying with skip.feature=%d, skip.cell=%d (original error: %s)",
      basename(mtx_path), skip.feature, skip.cell, msg
    ))
    
    Seurat::ReadMtx(mtx = mtx_path, cells = barcodes_path, features = genes_path,
                    feature.column = 1, skip.cell = skip.cell, skip.feature = skip.feature)
  })
}

parse_HTAPP <- function(sample_id, mtx_path, genes_path, barcodes_path, cell_index_df) {
  counts <- read_mtx_safe(mtx_path, genes_path, barcodes_path)
  
  # Read genes/barcodes with the same skip logic as whatever ReadMtx ended up using,
  # by re-deriving the offsets the same way (cheap: just re-check line counts)
  dims <- get_mtx_dims(mtx_path)
  n_genes_lines    <- length(readLines(genes_path))
  n_barcodes_lines <- length(readLines(barcodes_path))
  skip.feature <- max(n_genes_lines    - dims$n_features, 0)
  skip.cell    <- max(n_barcodes_lines - dims$n_barcodes, 0)
  
  genes    <- read.delim(genes_path,    header = FALSE, skip = skip.feature)$V1
  barcodes <- read.delim(barcodes_path, header = FALSE, skip = skip.cell)$V1
  
  stopifnot(length(genes) == nrow(counts), length(barcodes) == ncol(counts))
  
  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    rowData = data.frame(gene_id = genes, row.names = genes),
    colData = data.frame(barcode = barcodes, sample_id = sample_id,
                         row.names = barcodes)
  )
  
  sce <- filter_empty_droplets(sce)
  
  cd_base <- colData(sce) |> as.data.frame() |> rownames_to_column("old_cell_id")
  
  if (is.null(cell_index_df)) {
    cd <- cd_base |> mutate(cell_id = seq_len(n()))
  } else {
    cd <- cd_base |>
      inner_join(cell_index_df, by = c("old_cell_id" = "NAME", "sample_id")) |>
      dplyr::rename(cell_id = cell_index)
    if (nrow(cd) == 0) return(NULL)
  }
  
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}

#' Build the global NAME -> cell_index map from HTAPP L4 TSV files.
#' Each row is: sample_id | NAME | cell_index (row_number within sample)
build_htapp_cell_index <- function(l4_tsv_paths) {
  purrr::pmap_dfr(l4_tsv_paths, function(sample_id, processed_file) {
    read.table(processed_file, sep = "\t", header = TRUE,
               stringsAsFactors = FALSE,
               colClasses = "character")[-1, ] |>
      mutate(sample_id = sample_id, source_file = processed_file)  # keep file provenance
  }) |>
    group_by(sample_id) |>
    mutate(cell_index = row_number()) |>   # Continuous across all files for that sample, avoid duplicate index for a sample read from different files.
    ungroup() |>
    select(sample_id, NAME, cell_index)
}


# -----------------------------------------------------------------------------
# Parser: HTAN BU  (HTA3)
# Strategy: L3 dense CSV  +  L4 CSV cell filter
#
# L3: Dense gene × cell CSV (raw counts). One file per biospecimen.
#     Column names use "." instead of "-" in barcodes; rownames have version suffix.
# L4: CSV with SampleID and NAME columns covering multiple biospecimens.
#     NOTE: SampleID in L4 misses a leading digit vs Biospecimen — fix with regex.
# -----------------------------------------------------------------------------
parse_BU <- function(counts_path, biospecimen, cell_index_df) {
  raw <- read.csv(counts_path, row.names = 1, check.names = FALSE) |> as.matrix()
  
  # Restore 10X barcode format and strip gene version suffix
  colnames(raw) <- gsub("\\.([0-9]+)$", "-\\1", colnames(raw))
  rownames(raw) <- gsub("\\.\\d+$",     "",      rownames(raw))
  
  sce <- SingleCellExperiment(assays = list(counts = raw))
  sce$sample_id <- biospecimen
  
  cd <- colData(sce) |> as.data.frame() |> rownames_to_column("old_cell_id") |>
    inner_join(cell_index_df, by = c("old_cell_id" = "NAME", "sample_id"="SampleID")) |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}

#' Build cell index from BU L4 CSV files.
#' Fixes the known SampleID digit-drop bug (HTA3_8001_001 -> HTA3_8001_1001).
build_bu_cell_index <- function(l4_csv_paths) {
  purrr::map_dfr(l4_csv_paths, ~ read.csv(.x, colClasses = "character") |> slice(-1)) |>
    mutate(SampleID = sub("_(\\d+)$", "_1\\1", SampleID)) |>
    group_by(SampleID) |>
    mutate(cell_index = row_number()) |>
    ungroup() |>
    select(NAME, cell_index, SampleID)
}


# -----------------------------------------------------------------------------
# Parser: HTAN MSK  (HTA8)
# Strategy: L4 h5ad (contains raw counts + QC-filtered cells)
#
# WHY L4 AND NOT L3?
#   - MSK L4 h5ads cover 121 biospecimens; 3 have NO L3 file at all.
#   - The L4 merges all samples from a cohort: one read call vs 50-120 file reads.
#   - assay[1] is raw counts (integer), not normalised.
#   - The `patient` colData column encodes the internal MSK sample ID.
#     An external CSV maps internal IDs to HTAN biospecimen IDs.
#
# map_id_path: CSV with columns `sample` (MSK internal) and `sample_HTAN_ID`.
#              sample_HTAN_ID may be semicolon-separated (one MSK sample = several HTAN IDs).
# -----------------------------------------------------------------------------
parse_MSK <- function(h5ad_path, target_sample_ids, cell_index_df, map_id_path) {
  map_id <- read.csv(map_id_path)
  data   <- zellkonverter::readH5AD(h5ad_path, reader = "R", use_hdf5 = TRUE)
  
  if (!"patient" %in% colnames(as.data.frame(colData(data)))) return(NULL)
  
  # Keep only the `patient` column, then attach HTAN IDs
  colData(data) <- colData(data)[, "patient", drop = FALSE]
  data <- data |>
    tidyr::separate(".cell", c("sample", "barcode"), sep = "_(?=[^_]+$)", remove = FALSE) |>
    left_join(map_id, by = "sample") |>
    dplyr::rename(sample_id = sample_HTAN_ID)
  
  keep <- which(colData(data)$sample_id %in% target_sample_ids)
  if (length(keep) == 0) return(NULL)
  data <- data[, keep]
  
  # Restrict to first assay (raw counts), rename to "counts"
  assay_name <- names(assays(data))[1]
  assays(data) <- assays(data)[assay_name]
  assayNames(data) <- "counts"
  
  # Apply global cell index
  cd <- colData(data) |> as.data.frame() |> rownames_to_column("old_cell_id") |>
    inner_join(cell_index_df, by = c("sample_id", "old_cell_id" = ".cell"))
  
  if (nrow(cd) == 0) return(NULL)
  data <- data[, cd$old_cell_id]
  stopifnot(ncol(data) == nrow(cd))
  finalise_sce(data, setNames(as.character(cd$cell_index), cd$old_cell_id))
}


# -----------------------------------------------------------------------------
# Parser: HTAN OHSU  (HTA9)
# Strategy: L3 10X HDF5  (1:1 per biospecimen, no L4 available)
#
# Files: numeric IDs like 0000368073.h5 (CellRanger HDF5 format).
# All cells are used (no L4 cell filter exists in the breast metadata).
# cell_index_df: sample-specific NAME/cell_index table from ohsu_cell_index_df.
# -----------------------------------------------------------------------------
parse_OHSU <- function(h5_path, sample_id, cell_index_df) {
  sce <- DropletUtils::read10xCounts(h5_path, type = "HDF5", col.names = TRUE)
  sce$sample_id <- sample_id
  
  sce <- filter_empty_droplets(sce)
  
  cd <- colData(sce) |> as.data.frame() |> rownames_to_column("old_cell_id") |>
    inner_join(cell_index_df, by = c("old_cell_id" = "NAME", "sample_id")) |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}

#' Build cell index for OHSU samples by reading barcodes directly from HDF5.
#' Uses rhdf5 to avoid loading the full count matrix.
#' Returns a tibble with sample_id | NAME (barcode) | cell_index.
build_ohsu_cell_index <- function(metadata) {
  purrr::pmap_dfr(metadata, function(sample_id, full_path) {
    barcodes <- rhdf5::h5read(full_path, "/matrix/barcodes")
    tibble::tibble(
      sample_id  = sample_id,
      NAME       = barcodes,
      cell_index = seq_along(barcodes)
    )
  })
}


# -----------------------------------------------------------------------------
# Parser: HTAN WUSTL  (HTA12)
# Strategy: L4 per-sample Seurat RDS as sole data source (L3 MTX not used).
#
# L4: One .rds Seurat object per biospecimen (some have multiple RDS files —
#     pick the first after sorting by filename for reproducibility).
#     Active assay is SCT; RNA assay counts are extracted explicitly.
# cell_index_df: sample-specific NAME/cell_index table from wustl_cell_index_df.
# -----------------------------------------------------------------------------

#' Build cell index for WUSTL from per-sample L4 Seurat RDS files.
#' Barcodes are sorted within each sample for a deterministic cell_index.
#' Returns tibble: sample_id | NAME (barcode) | cell_index
build_wustl_cell_index <- function(metadata) {
  purrr::pmap_dfr(metadata, function(sample_id, l4_rds_path) {
    seurat_obj <- readRDS(l4_rds_path)
    barcodes   <- sort(Seurat::Cells(seurat_obj))
    tibble::tibble(
      sample_id  = sample_id,
      NAME       = barcodes,
      cell_index = seq_along(barcodes)
    )
  })
}

#' Parse a single WUSTL sample from its L4 Seurat RDS.
#' Extracts RNA assay counts (not SCT), assigns sample_id, then maps
#' cell barcodes to deterministic indices via cell_index_df.
parse_WUSTL <- function(sample_id, l4_rds_path, cell_index_df) {
  seurat_obj <- readRDS(l4_rds_path)
  counts     <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    rowData = data.frame(gene_id = rownames(counts), row.names = rownames(counts)),
    colData = data.frame(barcode   = colnames(counts),
                         sample_id = sample_id,
                         row.names = colnames(counts))
  )
  
  cd <- colData(sce) |> as.data.frame() |> tibble::rownames_to_column("old_cell_id") |>
    dplyr::inner_join(cell_index_df, by = c("old_cell_id" = "NAME"), "sample_id") |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}


# -----------------------------------------------------------------------------
# Parser: HTAN Duke  (HTA6)
# Strategy: L4 single merged RDS (22 samples) as primary source;
#           L3 MTX triplet as fallback for samples absent from L4.
#
# L4: ONE RDS Seurat object (HTAN_BIOSPECIMEN_ID column per cell).
#     Used directly as the data source — no L3 MTX loading for these samples.
#     Cell index = sorted barcodes within each sample (reproducible).
# L3: MTX triplet fallback for samples present in L3 but missing from L4.
#     Cell index = sequential over barcodes file (via duke_cell_index_df).
# -----------------------------------------------------------------------------

#' Build the full Duke cell index:
#'  - L4 Seurat samples : sort barcodes within sample → deterministic cell_index
#'  - L3-only samples   : read barcodes file → sequential cell_index
#' Returns tibble: sample_id | NAME (barcode) | cell_index
build_duke_cell_index <- function(l4_rds_path, l3_metadata) {
  seurat_obj <- readRDS(l4_rds_path)
  meta       <- seurat_obj@meta.data |> tibble::rownames_to_column(".cell")
  
  if (!"HTAN_BIOSPECIMEN_ID" %in% colnames(meta))
    stop("Expected HTAN_BIOSPECIMEN_ID column in Duke Seurat metadata. ",
         "Found: ", paste(colnames(meta), collapse = ", "))
  
  l4_index <- meta |>
    dplyr::rename(sample_id = HTAN_BIOSPECIMEN_ID, NAME = .cell) |>
    dplyr::select(sample_id, NAME) |>
    dplyr::group_by(sample_id) |>
    dplyr::arrange(NAME, .by_group = TRUE) |>
    dplyr::mutate(cell_index = dplyr::row_number()) |>
    dplyr::ungroup()
  
  l4_samples <- unique(l4_index$sample_id)
  missing    <- l3_metadata |>
    dplyr::filter(!sample_id %in% l4_samples) |>
    dplyr::select(sample_id, barcodes_path) |>
    dplyr::distinct()
  
  if (nrow(missing) == 0) return(l4_index)
  
  l3_index <- purrr::pmap_dfr(missing, function(sample_id, barcodes_path) {
    barcodes <- read.delim(barcodes_path, header = FALSE)$V1
    tibble::tibble(
      sample_id  = sample_id,
      NAME       = barcodes,
      cell_index = seq_along(barcodes)
    )
  })
  
  dplyr::bind_rows(l4_index, l3_index)
}

#' Parse all L4-covered Duke samples directly from the merged Seurat RDS.
#' Splits by HTAN_BIOSPECIMEN_ID, converts each subset to a lean SCE using
#' pre-built cell_index_df for deterministic column naming.
#' Returns a named list of SCEs (names = sample_id).
parse_Duke_seurat <- function(l4_rds_path, cell_index_df) {
  seurat_obj <- readRDS(l4_rds_path)
  meta       <- seurat_obj@meta.data |> tibble::rownames_to_column(".cell")
  l4_samples <- unique(meta[["HTAN_BIOSPECIMEN_ID"]])
  
  purrr::set_names(l4_samples) |>
    purrr::map(function(sid) {
      cells   <- meta |> dplyr::filter(HTAN_BIOSPECIMEN_ID == sid) |> dplyr::pull(.cell)
      sub_seu <- seurat_obj[, cells]
      counts  <- Seurat::GetAssayData(sub_seu, assay = "RNA", layer = "counts")
      
      sce <- SingleCellExperiment(
        assays  = list(counts = counts),
        rowData = data.frame(gene_id = rownames(counts), row.names = rownames(counts)),
        colData = data.frame(barcode   = colnames(counts),
                             sample_id = sid,
                             row.names = colnames(counts))
      )
      
      cid <- cell_index_df |>
        dplyr::filter(sample_id == sid) |>
        dplyr::select(NAME, cell_index, sample_id)
      
      cd <- colData(sce) |> as.data.frame() |> tibble::rownames_to_column("old_cell_id") |>
        dplyr::inner_join(cid, by = c("old_cell_id" = "NAME", "sample_id")) |>
        dplyr::rename(cell_id = cell_index)
      
      if (nrow(cd) == 0) return(NULL)
      sce_sub <- sce[, cd$old_cell_id]
      finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
    })
}

#' Parse a single L3-only Duke sample (no L4 Seurat coverage).
#' Uses cell_index_df for barcode filtering and index assignment.
parse_Duke <- function(sample_id, mtx_path, genes_path, barcodes_path,
                       cell_index_df) {
  counts   <- Seurat::ReadMtx(mtx = mtx_path, cells = barcodes_path,
                              features = genes_path, feature.column = 1)
  genes    <- read.delim(genes_path,    header = FALSE)$V1
  barcodes <- read.delim(barcodes_path, header = FALSE)$V1
  
  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    rowData = data.frame(gene_id = genes, row.names = genes),
    colData = data.frame(barcode = barcodes, sample_id = sample_id,
                         row.names = barcodes)
  )
  
  sce <- filter_empty_droplets(sce)
  
  cd <- colData(sce) |> as.data.frame() |> tibble::rownames_to_column("old_cell_id") |>
    dplyr::inner_join(cell_index_df, by = c("old_cell_id" = "NAME", "sample_id")) |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}


# -----------------------------------------------------------------------------
# Parser: HTAN DFCI  (HTA5) 
# Strategy: L3 10X HDF5 (multiple biospecimens pooled per lane, no L4)
#
# -----------------------------------------------------------------------------
parse_DFCI <- function(h5_path, sample_id, cell_index_df) {
  sce <- DropletUtils::read10xCounts(h5_path, type = "HDF5", col.names = TRUE)
  sce$sample_id <- sample_id
  
  sce <- filter_empty_droplets(sce)
  
  cd <- colData(sce) |> as.data.frame() |> rownames_to_column("old_cell_id") |>
    inner_join(cell_index_df, by = c("old_cell_id" = "NAME", "sample_id")) |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}

#' Build cell index for DFCI pooled samples from HDF5 barcodes.
#' Returns a tibble with sample_id | NAME (barcode) | cell_index.
build_dfci_cell_index <- function(metadata) {
  purrr::pmap_dfr(metadata, function(sample_id, full_path) {
    barcodes <- rhdf5::h5read(full_path, "/matrix/barcodes")
    tibble::tibble(
      sample_id    = sample_id,
      NAME       = barcodes,
      cell_index = seq_along(barcodes)
    )
  })
}

# -----------------------------------------------------------------------------
# Parser: HTAN CHOP  (HTA4)
# Strategy: L4 per-sample Seurat RDS when available (NBL cohort);
#           L3 MTX triplet fallback for samples with no L4 (ped-glioma cohort).
#
# L4: seurat_objects_count_based/seurat_object_NB_7767_*.rds
#     One .rds per biospecimen (single Biospecimen ID in metadata).
#     These are the raw-count Seurat objects submitted by CHOP for the NBL
#     cohort.  The other L4 RDS files (Global_seurat_object.rds, merged
#     patient/ped-glioma objects) are integration/QC artefacts — ignored here.
# L3: Standard CellRanger triplet per sample directory:
#       barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz
#     Files use plain basenames (no per-sample prefix unlike Duke/HTAPP).
#     Rows with Biospecimen == "0 Biospecimens" are backup artefacts — excluded.
#
# Cell index:
#   L4 samples : sorted barcodes from Seurat object → deterministic cell_index
#   L3 samples : sequential barcodes from barcodes.tsv.gz
# -----------------------------------------------------------------------------

#' Build per-sample CHOP cell index from L4 per-sample Seurat RDS.
#' Called once per sample (patterned target).
#' Returns tibble: sample_id | NAME (barcode) | cell_index
build_chop_l4_cell_index <- function(metadata) {
  purrr::pmap_dfr(metadata, function(sample_id, l4_rds_path) {
    seurat_obj <- readRDS(l4_rds_path)
    seurat_obj <- UpdateSeuratObject(seurat_obj)
    barcodes   <- sort(Seurat::Cells(seurat_obj))
    tibble::tibble(
      sample_id  = sample_id,
      NAME       = barcodes,
      cell_index = seq_along(barcodes)
    )
  })
}

#' Build CHOP cell index for L3-only samples (no L4 RDS) from barcodes files.
#' Returns tibble: sample_id | NAME (barcode) | cell_index
build_chop_l3_cell_index <- function(metadata) {
  purrr::pmap_dfr(metadata, function(sample_id, barcodes_path) {
    barcodes <- read.delim(barcodes_path, header = FALSE)$V1
    tibble::tibble(
      sample_id  = sample_id,
      NAME       = barcodes,
      cell_index = seq_along(barcodes)
    )
  })
}

#' Parse a single CHOP sample from its L4 per-sample Seurat RDS.
#' Extracts RNA assay counts, assigns sample_id, maps barcodes to the
#' deterministic indices supplied in cell_index_df.
#' No emptyDrops filter — L4 Seurat objects are already QC-filtered.
parse_CHOP_seurat <- function(sample_id, l4_rds_path, cell_index_df) {
  seurat_obj <- readRDS(l4_rds_path)
  seurat_obj <- UpdateSeuratObject(seurat_obj)
  counts     <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    rowData = data.frame(gene_id = rownames(counts), row.names = rownames(counts)),
    colData = data.frame(barcode   = colnames(counts),
                         sample_id = sample_id,
                         row.names = colnames(counts))
  )
  
  cd <- colData(sce) |> as.data.frame() |> tibble::rownames_to_column("old_cell_id") |>
    dplyr::inner_join(cell_index_df, by = c("old_cell_id" = "NAME")) |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}

#' Parse a single CHOP L3-only sample from a CellRanger MTX triplet.
#' Applies emptyDrops filtering before joining to cell_index_df.
parse_CHOP <- function(sample_id, mtx_path, genes_path, barcodes_path,
                       cell_index_df) {
  counts   <- Seurat::ReadMtx(mtx = mtx_path, cells = barcodes_path,
                              features = genes_path, feature.column = 1)
  genes    <- read.delim(genes_path,    header = FALSE)$V1
  barcodes <- read.delim(barcodes_path, header = FALSE)$V1
  
  genes          <- strip_gene_version_suffix(genes)
  rownames(counts) <- genes
  
  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    rowData = data.frame(gene_id = genes, row.names = genes),
    colData = data.frame(barcode = barcodes, sample_id = sample_id,
                         row.names = barcodes)
  )
  
  sce <- filter_empty_droplets(sce)
  
  cd <- colData(sce) |> as.data.frame() |> tibble::rownames_to_column("old_cell_id") |>
    dplyr::inner_join(cell_index_df, by = c("old_cell_id" = "NAME", "sample_id")) |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}

# -----------------------------------------------------------------------------
# Parser: HTAN Stanford  (HTA10)
# Strategy: L3 MTX triplet only (no L4 available)
#
# Files: {SamplePrefix}_barcodes.tsv.gz / _features.tsv.gz / _matrix.mtx.gz
# Some biospecimens have both a primary run and an "-R0" re-run that share
# the same Biospecimen ID.  The "-R0" files are excluded at the metadata-prep
# step (filter Filename_basename containing "R0") so each biospecimen maps
# to exactly one MTX triplet.
# -----------------------------------------------------------------------------

#' Build cell index for Stanford samples from barcodes files.
#' Returns tibble: sample_id | NAME (barcode) | cell_index
build_stanford_cell_index <- function(metadata) {
  purrr::pmap_dfr(metadata, function(sample_id, barcodes_path) {
    barcodes <- read.delim(barcodes_path, header = FALSE)$V1
    tibble::tibble(
      sample_id  = sample_id,
      NAME       = barcodes,
      cell_index = seq_along(barcodes)
    )
  })
}

#' Parse a single Stanford sample from a CellRanger MTX triplet.
#' Applies emptyDrops filtering; no L4 cell filter is available.
parse_stanford <- function(sample_id, mtx_path, genes_path, barcodes_path,
                           cell_index_df) {
  counts   <- Seurat::ReadMtx(mtx = mtx_path, cells = barcodes_path,
                              features = genes_path, feature.column = 1)
  genes    <- read.delim(genes_path,    header = FALSE)$V1
  barcodes <- read.delim(barcodes_path, header = FALSE)$V1
  
  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    rowData = data.frame(gene_id = genes, row.names = genes),
    colData = data.frame(barcode = barcodes, sample_id = sample_id,
                         row.names = barcodes)
  )
  
  sce <- filter_empty_droplets(sce)
  
  cd <- colData(sce) |> as.data.frame() |> tibble::rownames_to_column("old_cell_id") |>
    dplyr::inner_join(cell_index_df, by = c("old_cell_id" = "NAME", "sample_id")) |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}

# -----------------------------------------------------------------------------
# Parser: HTAN Vanderbilt  (HTA11)
# Strategy: L3 dense CSV, cell-type-partitioned or multi-batch cbind
#
# Two CSV sub-types share the same parser:
#   YX series:  {batch}-YX-{N}_epi.csv + {batch}-YX-{N}_immune.csv
#               Always 2 files per biospecimen (epithelial + immune compartments).
#               → cbind the two cell partitions.
#   as series:  {batch}as{N}.csv — 1–4 CSVs from different study batches
#               mapping to the same Biospecimen ID (different cells per batch).
#               → cbind all files after gene-set intersection.
#
# Format: gene × cell dense CSV (rows = gene symbols, cols = cell IDs).
# No L4 available; all called cells are retained (no emptyDrops needed).
# Barcode collisions across files are resolved with a ___f{i} suffix.
#
# MTX sub-group (P9142 only, 2 biospecimens): parsed with parse_stanford().
# Lau MTX group (~13 biospecimens): POOLED lanes — flag like DFCI.
# -----------------------------------------------------------------------------

#' Build cell index for Vanderbilt CSV samples.
#' metadata must have columns sample_id and csv_paths (a list-column of sorted paths).
#' Reads only column headers to get cell IDs without loading full matrices.
#' Returns tibble: sample_id | NAME (cell ID) | cell_index
build_vanderbilt_csv_cell_index <- function(metadata) {
  purrr::pmap_dfr(metadata, function(sample_id, csv_paths) {
    # cell_ids_per_file <- lapply(csv_paths, function(p) {
    #   rownames(read.csv(resolve_csv_path(p), row.names = 1, nrows = 0, check.names = FALSE))
    # })
    
    cell_ids_per_file <- lapply(csv_paths, function(p) {
      # Format: cell × gene (rows = cells, cols = genes).
      # Cell IDs are in the first column of each data row (not the header).
      # Read lines directly and extract the first comma-delimited field,
      # skipping the header row. Avoids loading the full count matrix.
      con <- file(resolve_csv_path(p), open = "r")
      on.exit(close(con))
      lines <- readLines(con)
      vapply(lines[-1L], function(l) {
        end <- regexpr(",", l, fixed = TRUE)
        if (end == -1L) l else substr(l, 1L, end - 1L)
      }, character(1L), USE.NAMES = FALSE)
    })
    
    all_ids <- unlist(cell_ids_per_file)
    n_dup   <- sum(duplicated(all_ids))
    if (n_dup > 0) {
      message("build_vanderbilt_csv_cell_index [", sample_id, "]: dropping ",
              n_dup, " duplicated cell ID(s) (keeping first occurrence).")
      all_ids <- all_ids[!duplicated(all_ids)]
    }
    tibble::tibble(
      sample_id  = sample_id,
      NAME       = all_ids,
      cell_index = seq_along(all_ids)
    )
  })
}

#' Parse a single Vanderbilt sample from 1–N dense CSVs.
#' Handles single CSVs, epi+immune pairs, and multi-batch as-series equally.
#' CSVs are assumed gene × cell (rows = gene symbols, cols = cell IDs).
#' If files are transposed (cell × gene), add t() in the read step.
parse_vanderbilt_csv <- function(sample_id, csv_paths, cell_index_df) {
  # mats <- lapply(csv_paths, function(p) {
  #   as.matrix(t(read.csv(resolve_csv_path(p), row.names = 1, check.names = FALSE)))
  # })
  
  mats <- lapply(csv_paths, function(p) {
    # Format: cell × gene (rows = cells, cols = genes).
    # Read without row.names = 1 to avoid errors on duplicate gene names.
    # First column = cell IDs; remaining columns = gene expression values.
    # Transpose to gene × cell after extracting cell IDs and gene names.
    df    <- read.csv(resolve_csv_path(p), check.names = FALSE, row.names = NULL)
    cells <- df[[1L]]
    genes <- colnames(df)[-1L]
    m     <- t(as.matrix(df[, -1L, drop = FALSE]))  # → gene × cell
    colnames(m) <- cells
    # Deduplicate gene rows (keep first occurrence of each symbol)
    keep_genes <- !duplicated(genes)
    if (!all(keep_genes))
      message("parse_vanderbilt_csv [", sample_id, "]: dropping ",
              sum(!keep_genes), " duplicate gene(s) in ", basename(p), ".")
    m           <- m[keep_genes, , drop = FALSE]
    rownames(m) <- genes[keep_genes]
    m
  })
  
  # Intersect gene sets (should be identical within a series, but be defensive)
  common_genes <- Reduce(intersect, lapply(mats, rownames))
  if (length(common_genes) == 0)
    stop("parse_vanderbilt_csv: no common genes across CSV files for ", sample_id)
  if (length(common_genes) < nrow(mats[[1]]))
    warning("parse_vanderbilt_csv [", sample_id, "]: restricting to ",
            length(common_genes), " common genes (from ", nrow(mats[[1]]), ").")
  mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
  
  counts <- do.call(cbind, mats)
  
  # Drop duplicate cell IDs (keep first occurrence, matching cell index builder)
  dup_cells <- duplicated(colnames(counts))
  if (any(dup_cells)) {
    message("parse_vanderbilt_csv [", sample_id, "]: dropping ",
            sum(dup_cells), " duplicated cell ID(s) (keeping first occurrence).")
    counts <- counts[, !dup_cells, drop = FALSE]
  }
  
  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    rowData = data.frame(gene_id = common_genes, row.names = common_genes),
    colData = data.frame(sample_id = rep(sample_id, ncol(counts)), row.names = colnames(counts))
  )
  
  cd <- colData(sce) |> as.data.frame() |> tibble::rownames_to_column("old_cell_id") |>
    dplyr::inner_join(cell_index_df, by = c("old_cell_id" = "NAME", "sample_id")) |>
    dplyr::rename(cell_id = cell_index)
  
  if (nrow(cd) == 0) return(NULL)
  sce_sub <- sce[, cd$old_cell_id]
  stopifnot(ncol(sce_sub) == nrow(cd))
  finalise_sce(sce_sub, setNames(cd$cell_id, cd$old_cell_id))
}

#' Parse a single Vanderbilt MTX sample (P9142 sub-group).
#' Delegates to parse_stanford — same CellRanger triplet format, same
#' emptyDrops filter, no L4. Named separately for registry clarity.
parse_vanderbilt <- function(...) parse_stanford(...)

# -----------------------------------------------------------------------------
# Center registry
# -----------------------------------------------------------------------------
#
# Each entry describes HOW to prepare the metadata before calling the parser,
# and WHICH parser to call. Use this for documentation and dispatch logic in
# the targets scripts below.
#
# strategy:       short label for the parsing approach
# l3_file_types:  regex patterns that identify the relevant L3 files
# l4_file_types:  regex patterns that identify the relevant L4 files (or NULL)
# parse_fn:       the function to call
# note:           any known gotchas for this center
# -----------------------------------------------------------------------------
CENTRE_REGISTRY <- list(
  
  "HTAN HTAPP" = list(
    strategy      = "L3_MTX + L4_TSV_filter",
    l3_file_types = list(
      mtx_path      = "_matrix\\.mtx\\.gz$",
      genes_path    = "_features\\.tsv\\.gz$",
      barcodes_path = "_barcodes\\.tsv\\.gz$"
    ),
    l4_file_types = list(processed_file = "_L4\\.tsv$"),
    l4_level      = "Level 4",
    parse_fn      = parse_HTAPP,
    channel_split = TRUE,   # biospecimens with >4 files are split per channel
    note          = paste(
      "Multi-channel biospecimens: sample_id = paste(Biospecimen, channel, sep='___').",
      "L4 TSV is a batch file covering many biospecimens; build global cell index first.",
      "Breast has 'Unknown' h5ads in data_integration_MBC/ — IGNORE (integration outputs,",
      "strict subset of MTX triplets)."
    )
  ),
  
  "HTAN BU" = list(
    strategy      = "L3_CSV + L4_CSV_filter",
    l3_file_types = list(counts_path = "raw.*\\.csv$"),
    l4_file_types = list(processed_file = "\\.csv$"),
    l4_level      = "Level 4",
    parse_fn      = parse_BU,
    channel_split = FALSE,
    note          = paste(
      "L4 CSV covers all biospecimens in one file (filter by SampleID column).",
      "Known bug: SampleID misses a leading digit vs Biospecimen — fix with",
      "sub('_(\\d+)$', '_1\\1', SampleID).",
      "Gene names are symbols — convert to Ensembl with convert_gene_to_ensembl()."
    )
  ),
  
  "HTAN MSK" = list(
    strategy      = "L4_H5AD (raw counts embedded)",
    l3_file_types = NULL,  # L3 not needed
    l4_file_types = list(h5ad_path = "\\.h5ad$"),
    l4_level      = "Level 4",
    parse_fn      = parse_MSK,
    channel_split = FALSE,
    note          = paste(
      "L4 h5ads cover 121 biospecimens; 3 have no L3 file — L4 is the ONLY source.",
      "Three separate cohorts in Lung (Laughney/Bhatt, Chan, Glasner) each have their",
      "own h5ad. Requires adata_sample_id_htan_id_map.csv to map MSK internal IDs",
      "to HTAN biospecimen IDs. assay[1] is raw counts."
    )
  ),
  
  "HTAN OHSU" = list(
    strategy      = "L3_10X_HDF5 (1:1 per biospecimen)",
    l3_file_types = list(h5_path = "\\.h5$"),
    l4_file_types = NULL,
    l4_level      = NULL,
    parse_fn      = parse_OHSU,
    channel_split = FALSE,
    note          = "No L4 in Breast metadata. All cells used. File naming: 0000368073.h5."
  ),
  
  "HTAN WUSTL" = list(
    strategy      = "L3_MTX + L4_RDS_per_sample",
    l3_file_types = list(
      mtx_path      = "-matrix\\.mtx\\.gz$",
      genes_path    = "-features\\.tsv\\.gz$",
      barcodes_path = "-barcodes\\.tsv\\.gz$"
    ),
    l4_file_types = list(l4_rds_path = "\\.rds$"),
    l4_level      = "Level 4",
    parse_fn      = parse_WUSTL,
    channel_split = FALSE,
    note          = paste(
      "Files are in scsnRNA-Seq_level_3_atac_tumor/ — multiome experiment.",
      "features.tsv.gz contains only RNA features; MTX read is RNA only.",
      "File basename (HT###B1-S1H#) ≠ HTAN biospecimen ID (HTA12_xxx_1);",
      "use Biospecimen column in files_metadata to build the mapping."
    )
  ),
  
  "HTAN Duke" = list(
    strategy      = "L3_MTX + L4_single_merged_RDS",
    l3_file_types = list(
      mtx_path      = "-matrix\\.mtx\\.gz$",
      genes_path    = "-features\\.tsv\\.gz$",
      barcodes_path = "-barcodes\\.tsv\\.gz$"
    ),
    l4_file_types = list(l4_rds_path = "\\.rds$"),
    l4_level      = "Level 4",
    parse_fn      = parse_Duke,
    channel_split = FALSE,
    note          = paste(
      "L4 is ONE Seurat RDS covering all 22 Breast biospecimens.",
      "Extract per-sample cell barcodes using orig.ident == biospecimen_id.",
      "Use build_duke_cell_filter() to build a named list before parallel processing."
    )
  ),
  
  "HTAN DFCI" = list(
    strategy      = "L3_10X_HDF5_POOLED (no demultiplexing available)",
    l3_file_types = list(h5_path = "_raw_feature_bc_matrix\\.h5$"),
    l4_file_types = NULL,
    l4_level      = NULL,
    parse_fn      = parse_DFCI,
    channel_split = FALSE,
    note          = paste(
      "POOLED: 3-4 HTAN biospecimen IDs share one 10X lane.",
      "No per-cell sample assignment available in HTAN metadata.",
      "Reads return pool_id (file basename) not a per-sample ID.",
      "Resolve with HTO demultiplexing or SNP-based tools (vireo/demuxlet)."
    )
  ),
  
  "HTAN CHOP" = list(
    strategy      = "L4_RDS_per_sample (NBL) + L3_MTX_fallback (ped-glioma)",
    l3_file_types = list(
      mtx_path      = "^matrix\\.mtx\\.gz$",
      genes_path    = "^features\\.tsv\\.gz$",
      barcodes_path = "^barcodes\\.tsv\\.gz$"
    ),
    l4_file_types = list(l4_rds_path = "seurat_objects_count_based/.*\\.rds$"),
    l4_level      = "Level 4",
    parse_fn      = parse_CHOP,
    channel_split = FALSE,
    note          = paste(
      "L4 per-sample RDS only from seurat_objects_count_based/ (NBL cohort).",
      "Other L4 RDS files (Global_seurat_object, merged patient/ped-glioma) are",
      "integration artefacts — ignored. L3 MTX triplet uses plain basenames",
      "(barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz); match by exact",
      "basename, not prefix. Rows with Biospecimen == '0 Biospecimens' are",
      "_bak duplicates — filter before building metadata."
    )
  ),
  
  "HTAN Stanford" = list(
    strategy      = "L3_MTX (no L4 available)",
    l3_file_types = list(
      mtx_path      = "_matrix\\.mtx\\.gz$",
      genes_path    = "_features\\.tsv\\.gz$",
      barcodes_path = "_barcodes\\.tsv\\.gz$"
    ),
    l4_file_types = NULL,
    l4_level      = NULL,
    parse_fn      = parse_stanford,
    channel_split = FALSE,
    note          = paste(
      "Colon NOS. L3 MTX triplet only; no L4 available.",
      "Some biospecimens have a duplicate '-R0' re-run — exclude",
      "Filename_basename values containing 'R0' before building metadata",
      "to ensure each biospecimen maps to exactly one MTX triplet."
    )
  ),
  
  "HTAN Vanderbilt" = list(
    strategy      = "L3_CSV_partitioned (epi+immune cbind, or multi-batch cbind)",
    l3_file_types = list(csv_path = "\\.csv$"),
    l4_file_types = NULL,
    l4_level      = NULL,
    parse_fn      = parse_vanderbilt_csv,
    channel_split = FALSE,
    note          = paste(
      "Colon/Not Reported. Two CSV sub-types:",
      "YX series: always 2 files per biospecimen (_epi + _immune) — cbind cell compartments.",
      "as series: 1–4 CSVs from different study batches per biospecimen — cbind all.",
      "P9142 MTX (2 biospecimens, single IDs) parsed with parse_stanford().",
      "Lau MTX group (~13 biospecimens, comma-separated IDs) is POOLED — flag like DFCI.",
      "CSVs are gene × cell dense matrices; gene names are symbols → convert_gene_to_ensembl()."
    )
  )
)

#' Print a human-readable summary of the registry
summarise_registry <- function() {
  purrr::iwalk(CENTRE_REGISTRY, function(entry, center) {
    cat(sprintf("\n%-20s | %-40s\n", center, entry$strategy))
    cat(sprintf("%-20s | %s\n", "", entry$note))
  })
}
