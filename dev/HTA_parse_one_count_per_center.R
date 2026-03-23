library(dplyr)
library(tidySingleCellExperiment)
library(SummarizedExperiment)
library(stringr)
library(tidyr)
library(purrr)
library(zellkonverter)
library(tibble)

file_path = "/vast/scratch/users/shen.m/synapse_data/lung/counts/"

file_metadata <- read.csv("/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv",
                          sep = "\t", na.strings = c("NA",""), header = TRUE) |> as_tibble() |> 
  # THIS IS ASSIGNED WRONG TO BIOSPECIMEN. HTAN PHASE1 IS SOOOO COMPLEX!
  dplyr::filter(Filename != "single_cell_RNAseq_level_4_lung/lung_HTA1_203_332102_ch1_L4.tsv") |>
  mutate(
    Filename_basename = basename(Filename),
    full_path = file.path(file_path, Filename_basename))

# =============================================================================
# Atlas MSK (HTA8)
# =============================================================================
map_id <- read.csv("/vast/scratch/users/shen.m/synapse_data/lung/counts/adata_sample_id_htan_id_map.csv", header= T)

HTA8 <- file_metadata |> filter(Atlas.Name=="HTAN MSK") |> filter(Level == "Level 4", File.Format == "hdf5") |> 
   mutate(sce = map(full_path, ~{
    data = zellkonverter::readH5AD(.x, reader= "R", use_hdf5 = T)
    
    if (!"patient" %in% colnames(as.data.frame(colData(data)))) return(NULL)
    
    # Only keep original counts
    assay_name <- assays(data) |> names() |> magrittr::extract2(1)
    assays(data) <- assays(data)[assay_name]
    
    # Rename to counts
    assayNames(data)[assayNames(data) == assay_name] <- "counts"
    
    data = data |> tidyr::separate(".cell", c("sample", "barcode"), sep = "_(?=[^_]+$)",  remove=FALSE) |> 
      left_join(map_id, by = "sample", copy=TRUE) |> 
      dplyr::rename(sample_id = sample_HTAN_ID)
    
    
  }, .progress = T)) 

# Create a sample_level metadata from list of existing h5ads
h5ad_list <- list.files("/vast/scratch/users/shen.m/synapse_data/lung/counts/", pattern = ".h5ad", full.name = TRUE) |> 
  purrr::map(~ {
    data <- zellkonverter::readH5AD(.x, reader= "R", use_hdf5 = T)
    data$full_path <- .x
    data
  }, .progress = T
    )

h5ad_metadata <- h5ad_list |> map( ~ colData(.x) |> as.data.frame() |> rownames_to_column(".cell"), .progress = T) |> 
  bind_rows() |> as_tibble() |> tidyr::separate(".cell", c("sample", "barcode"), sep = "_(?=[^_]+$)",  remove=FALSE) |>
  left_join(map_id, by = "sample", copy=TRUE) |> 
  dplyr::rename(sample_id = sample_HTAN_ID) |> 
  filter(!is.na(sample_id)) |> 
  separate_rows(sample_id, sep = ";") |>
  distinct(sample_id, full_path) |> 
  add_count(sample_id)

h5ad_metadata

# SAVE TO DISK PIPELINE IS DONE SEPARATELY USING TARGETS
pipeline <- souce("/home/users/allstaff/shen.m/git_control/cellNexus/dev/HTAN_parse_MSK_center_split_bind_sce.R")

# =============================================================================
# Atlas HTAPP (HTA1) 
# =============================================================================
parse_HTAPP <- function(sample, mtx_path, genes_path, barcodes_path, processed_path) {
  counts <- ReadMtx(mtx = mtx_path, 
                    cells = barcodes_path, 
                    features = genes_path, 
                    # Use Ensemble 
                    feature.column = 1)
  
  genes <- read.delim(genes_path, header = FALSE)$V1 # Extract Ensemble ID 
  barcodes <- read.delim(barcodes_path, header = FALSE)$V1 
  filtered_df <- read.table(processed_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)[-1,]
  filtered_cell_id <- filtered_df |> pull(NAME)
  
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    rowData = data.frame(gene_id = genes) |> tibble::column_to_rownames("gene_id"),
    colData = data.frame(barcode = barcodes,
                         #counts_path = basename(mtx_path),
                         sample_id = sample)
    
  )
  
  sce_subset <- sce[, colnames(sce) %in% filtered_cell_id]
  
  colnames(sce_subset) <- paste0(sce_subset$barcode, "___", sample)
  
  sce_subset
  
  
  
}

HTAPP = file_metadata |> filter(Atlas.Name=="HTAN HTAPP") |> 
  #filter(Biospecimen == "HTA1_274_4891101")|>
  arrange(desc(Biospecimen)) |> 
  head(8) |> 
  mutate(
    file_type = case_when(
      str_detect(Filename_basename, "_matrix.mtx.gz$") ~ "mtx_path",
      str_detect(Filename_basename, "_features.tsv.gz$") ~ "genes_path",
      str_detect(Filename_basename, "_barcodes.tsv.gz$") ~ "barcodes_path",
      str_detect(Filename_basename, "_L4.tsv$") ~ "processed_file"
    ),
    channel_number = Filename_basename |>
      str_replace_all("ch(?=[0-9]+)", "channel") |>
      str_extract("channel[0-9]+")
  ) |> 
  group_split(Biospecimen) |>
  map_dfr(~ {
    subgroup <- .x
    if (nrow(subgroup) > 4) {
      subgroup <- subgroup |>
        mutate(sample_id = paste(Biospecimen, channel_number, sep = "___"))
    } else {
      subgroup <- subgroup |>
        mutate(sample_id = Biospecimen)
    }
    subgroup
  }) |>
  dplyr::select(full_path, Biospecimen, sample_id, file_type) |>
  pivot_wider(
    id_cols = c(Biospecimen, sample_id),
    names_from = file_type,
    values_from = full_path
  ) |> 
  mutate(sce = pmap(
    list(sample_id, mtx_path, genes_path, barcodes_path, processed_file), 
    ~ parse_HTAPP(..1, ..2, ..3, ..4, ..5), .progress = T))


HTAPP |> pull(sce)
HTAPP
map2(HTAPP$sample_id, HTAPP$sce, ~writeH5AD(.y, file = paste0("~/scratch/htan/hta/09-11-2025/counts/", .x, ".h5ad"), 
                                        compression = "gzip"), .progress = T)

# =============================================================================
# Atlas BU (HTA3) 
# =============================================================================

parse_BU <- function(counts_path) {
  counts <- read.csv(counts_path, )
}

# This implies that Level 4 data contains every Biospecimen
BU_processed_df <- file_metadata |> filter(Atlas.Name=="HTAN BU") |>filter( Biospecimen |> str_detect(",") ) |> 
  pull(full_path) |> map_dfr(~ read.csv(.x) |> dplyr::slice(-1))

BU <- file_metadata |> filter(Atlas.Name=="HTAN BU") |> 
  filter(full_path |> str_detect("raw") ) |> 
  head(2) |> 
  mutate(sce = map2(full_path, Biospecimen, ~ {
    # read counts
    raw_counts <- read.csv(.x, row.names = 1, check.names = FALSE) |> as.matrix()
    
    # fix names
    colnames(raw_counts) <- gsub("\\.([0-9]+)$", "-\\1", colnames(raw_counts)) # restore 10x barcode format
    rownames(raw_counts) <- gsub("\\.\\d+$", "", rownames(raw_counts))         # remove gene version suffix
    
    # create SCE
    sce = SingleCellExperiment(assays = list(counts = raw_counts))
    sce$sample_id <- .y
    
    
    filtered_cell_id = BU_processed_df |> pull(NAME)
    
    sce_subset <- sce[, colnames(sce) %in% filtered_cell_id]
    
    sce_subset
  }, .progress = T))

map2(BU$Biospecimen, BU$sce, ~writeH5AD(.y, file = paste0("~/scratch/htan/hta/09-11-2025/counts/", .x, ".h5ad"), 
                                        compression = "gzip"), .progress = T)


