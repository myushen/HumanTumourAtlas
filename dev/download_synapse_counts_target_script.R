library(targets)
library(tarchetypes)
library(crew.cluster)
library(synapser)
library(dplyr)
library(purrr)

# Lung tissue metadata is download from 
# https://data.humantumoratlas.org/explore?selectedFilters=%5B%7B%22value%22%3A%22scRNA-seq%22%2C%22group%22%3A%22assayName%22%2C%22count%22%3A1832%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22Synapse%22%2C%22group%22%3A%22downloadSource%22%2C%22count%22%3A656%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22Lung%22%2C%22group%22%3A%22organType%22%2C%22count%22%3A942%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22csv%22%2C%22group%22%3A%22FileFormat%22%2C%22count%22%3A606%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22hdf5%22%2C%22group%22%3A%22FileFormat%22%2C%22count%22%3A13%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22mtx%22%2C%22group%22%3A%22FileFormat%22%2C%22count%22%3A144%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22tsv%22%2C%22group%22%3A%22FileFormat%22%2C%22count%22%3A103%2C%22isSelected%22%3Afalse%7D%5D

# To do: what's streamline way of crawling metadata from HTAN API?

sample_meta <- read.csv("inst/extdata/samples_metadata_2025_10_21.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)
donor_meta <- read.csv("inst/extdata/donors_metadata_2025_10_21.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)
file_meta <- read.csv("inst/extdata/files_metadata_2025_10_21.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)

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

lung_counts_download_target_store = "/vast/scratch/users/shen.m/download_lung_counts_synapse_data_target_store"

# Define target for each synapse ID
tar_script({
  library(targets)
  library(tarchetypes)
  library(crew.cluster)
  library(synapser)
  library(dplyr)
  library(purrr)
  
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    cue = tar_cue(mode = "thorough"),
    
    workspace_on_error = TRUE,
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
        )))
  
  # Read file and get id
  get_download_ids <- function(df_path) {
    df = read.csv(df_path, sep = "\t", na.strings = c("NA",""), header = TRUE)
    ids = df |> pull(Synapse.Id) |> unique()
  }
  
  # Function to download data from Synapse
  download_synapse_data <- function(id, save_directory) {
    set.seed(123)
    
    if (!dir.exists(save_directory)) {
      dir.create(save_directory, recursive = TRUE)
    }
    
    synLogin(authToken =  Sys.getenv("SYNAPSE_TOKEN"))
    synGet(entity = id, downloadLocation = save_directory)
    print("saved successfully.. ")
  }
  
  list(
    tar_target(
      file_metadata,
      "/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/files_metadata_2025_10_21.tsv",
      deployment = "main"
    ),
    
    tar_target(
      synapse_id,
      get_download_ids(file_metadata),
      deployment = "main"
      ),
    tar_target(
      download_data,
      download_synapse_data(synapse_id,
                            "/vast/scratch/users/shen.m/synapse_data/lung/counts"),
      pattern = map(synapse_id)
    )
  )
        
  
}, script = paste0(lung_counts_download_target_store, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(lung_counts_download_target_store, "_target_script.R"), 
    store = lung_counts_download_target_store, 
    reporter = "summary"
  )
  
})

