library(targets)
library(tarchetypes)
library(crew.cluster)
library(synapser)
library(dplyr)
library(purrr)

# breast tissue metadata is download from 
# https://data.humantumoratlas.org/explore?selectedFilters=%5B%7B%22value%22%3A%22scRNA-seq%22%2C%22group%22%3A%22assayName%22%2C%22count%22%3A1832%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22Synapse%22%2C%22group%22%3A%22downloadSource%22%2C%22count%22%3A656%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22breast%22%2C%22group%22%3A%22organType%22%2C%22count%22%3A942%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22csv%22%2C%22group%22%3A%22FileFormat%22%2C%22count%22%3A606%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22hdf5%22%2C%22group%22%3A%22FileFormat%22%2C%22count%22%3A13%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22mtx%22%2C%22group%22%3A%22FileFormat%22%2C%22count%22%3A144%2C%22isSelected%22%3Afalse%7D%2C%7B%22value%22%3A%22tsv%22%2C%22group%22%3A%22FileFormat%22%2C%22count%22%3A103%2C%22isSelected%22%3Afalse%7D%5D

# To do: what's streamline way of crawling metadata from HTAN API?

sample_meta <- read.csv("inst/extdata/breast/samples_metadata.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)
donor_meta <- read.csv("inst/extdata/breast/donors_metadata.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)
file_meta <- read.csv("inst/extdata/breast/files_metadata.tsv", sep = "\t", na.strings = c("NA",""), header = TRUE) |> head(2)

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

breast_counts_download_target_store = "/vast/scratch/users/shen.m/download_breast_counts_synapse_data_target_store"

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
      "/home/users/allstaff/shen.m/git_control/HumanTumourAtlas/inst/extdata/breast/files_metadata.tsv",
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
                            "/vast/scratch/users/shen.m/synapse_data/breast/counts"),
      pattern = map(synapse_id),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      ) 
    )
  )
  
  
}, script = paste0(breast_counts_download_target_store, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(breast_counts_download_target_store, "_target_script.R"), 
    store = breast_counts_download_target_store, 
    reporter = "summary"
  )
  
})

