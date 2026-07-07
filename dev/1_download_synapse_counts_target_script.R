library(targets)
library(tarchetypes)
library(crew.cluster)
library(synapser)
library(dplyr)
library(purrr)

# Metadata is download from HTAN website
# https://humantumoratlas.org/explore?selectedFilters=%5B%7B%22value%22%3A%22scRNA-seq%22%2C%22group%22%3A%22assayName%22%2C%22count%22%3A126%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22Level+3%22%2C%22group%22%3A%22level%22%2C%22count%22%3A13%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22Level+4%22%2C%22group%22%3A%22level%22%2C%22count%22%3A48%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22Synapse%22%2C%22group%22%3A%22downloadSource%22%2C%22count%22%3A4562%2C%22isSelected%22%3Afalse%7D%5D
# Filters:
#    - scRNA-seq Assay
#    - Level 3 or 4
#    - Synapse Data Access
# Downloaded biospecimens (sample_meta), cases (donor_meta), and files (file_meta) metadata separately.

sample_meta <- read.csv("inst/extdata/samples_metadata_scRNAseq_synapse_level3_4.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)
donor_meta <- read.csv("inst/extdata/donors_metadata_scRNAseq_synapse_level3_4.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)
file_meta <- read.csv("inst/extdata/files_metadata_scRNAseq_synapse_level3_4.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)

# Log in to Synapse, create your own token, and save to Renviron
if (!nzchar(Sys.getenv("SYNAPSE_TOKEN", unset = ""))) {
  cat(
    'SYNAPSE_TOKEN="your_own_token"\n',
    file = path.expand("~/.Renviron"),
    append = TRUE
  )
}
# Test login works
# synLogin(authToken = Sys.getenv("SYNAPSE_TOKEN"))

# Start downloading using targets

download_htan_scRNAseq_l3_l4_synapse_target_store = "/vast/scratch/users/shen.m/download_htan_scRNAseq_l3_l4_synapse_target_store"

# Define target for each synapse ID
tar_script({
  library(targets)
  library(tarchetypes)
  library(crew.cluster)
  library(synapser)
  library(dplyr)
  library(purrr)
  library(crew)
  
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
    format = "qs",
    #debug = "dataset_id_sct_ea377f6e2d0ae2b7",
    workspace_on_error = TRUE,
    controller = controllers, 
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    ) 
  )
  
  # Read file and get id
  get_download_ids <- function(df_path) {
    df = read.csv(df_path, sep = "\t", na.strings = c("NA",""), header = TRUE)
    df |> distinct(Synapse.Id, Atlas.Name) |>
      mutate(Atlas.Name = sub(" ", "_", Atlas.Name))
  }
  
  # Function to download data from Synapse
  download_synapse_data <- function(id, center_id, save_directory) {
    set.seed(123)
    
    if (!dir.exists(save_directory)) {
      dir.create(save_directory, recursive = TRUE)
    }
    
    out_dir <- file.path(save_directory, center_id)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    synLogin(authToken = Sys.getenv("SYNAPSE_TOKEN"))
    synGet(entity = id, downloadLocation = out_dir)
    print("saved successfully.. ")
  }
  
  list(
    tar_target(
      file_metadata,
      "/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/files_metadata_scRNAseq_synapse_level3_4.tsv",
      deployment = "main"
    ),
    
    tar_target(
      synapse_df,
      get_download_ids(file_metadata),
      # |>
      #   filter(Atlas.Name == "HTAN_BU"),
      deployment = "main"
      ),
    
    tar_target(
      synapse_id,
      synapse_df$Synapse.Id,
      deployment = "main"
    ),
    
    tar_target(
      atlas_name,
      synapse_df$Atlas.Name,
      deployment = "main"
    ),
    
    tar_target(
      download_data,
      download_synapse_data(synapse_id,
                            atlas_name,
                            "/vast/scratch/users/shen.m/synapse_data"),
      pattern = map(synapse_id, atlas_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      ) 
    )
  )
        
  
}, script = paste0(download_htan_scRNAseq_l3_l4_synapse_target_store, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(download_htan_scRNAseq_l3_l4_synapse_target_store, "_target_script.R"), 
    store = download_htan_scRNAseq_l3_l4_synapse_target_store, 
    reporter = "summary"
  )
  
})

tar_progress_branches(store = download_htan_scRNAseq_l3_l4_synapse_target_store)


