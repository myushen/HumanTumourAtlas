library(cellNexus)
library(dplyr)
library(magrittr)
library(purrr)
library(zellkonverter)
library(tibble)
library(stringr)
x = get_metadata(cloud_metadata = get_metadata_url("hta_metadata.0.1.0.parquet"))

# YOU WOULD EXPECT HTA ATLAS HERE: 
x |> dplyr::count(atlas_id) 

# TEST ONE SAMPLE
sce = x |> get_single_cell_experiment()

dir = "/vast/scratch/users/shen.m/r_cache/R/cellNexus/hta/09-11-2025//counts/"

files = list.files(dir, pattern = ".h5ad", full.names = T)

tbl = map(files, ~readH5AD(.x, reader = "R", use_hdf5 = T) |> rownames() |> head(1), .progress = T) |> 
  set_names(files) |> enframe() |> mutate(value = unlist(value))

tbl |> filter(!value |> str_detect("ENSG"))