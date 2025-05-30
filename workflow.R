lib_directory <- paste0(workflow_directory, "lib/")

# Pre-process stage
source(paste0(lib_directory, "clean_counts.R"))
source(paste0(lib_directory, "clean_clusters.R"))

# Search stage
source(paste0(lib_directory, "batch_cluster_sequences.R"))

# Identify stage
source(paste0(lib_directory, "consolidate_blast_results.R"))
source(paste0(lib_directory, "top_blast_results_by_cluster.R"))

# Annotation stage
source(paste0(lib_directory, "normalize_name.R"))
source(paste0(lib_directory, "get_gbif_accepted_name.R"))
source(paste0(lib_directory, "apply_top_match_to_clusters.R"))
source(paste0(lib_directory, "get_weather_data.R"))
source(paste0(lib_directory, "apply_match_to_counts.R"))

# For GBIF DNA-derived dataset
source(paste0(lib_directory, "create_sample_cluster_dataset.R"))

# Taxonomy and Traits
source(paste0(lib_directory, "create_organism_dataset.R"))
source(paste0(lib_directory, "get_fungal_traits.R"))

# Filter for Fungi
source(paste0(lib_directory, "filter_for_fungi.R"))
source(paste0(lib_directory, "filter_for_wood_saprotroph.R"))



run_workflow <- function(batch_name, scata_dataset_name, scata_job_id, prompt=TRUE) {
  print(paste("[MycoPins Workflow] Batch:", batch_name,
              "; SCATA dataset:", scata_dataset_name,
              "; SCATA job:", scata_job_id))

  configuration <- mycopins_config(batch_name, scata_dataset_name, scata_job_id)

  mycopins_preprocess(configuration)
  if (prompt) {
    readline(prompt="Press [Enter] to continue")
  }

  mycopins_search(configuration)
  if (prompt) {
    readline(prompt="Press [Enter] to continue")
  }

  mycopins_identify(configuration)
  mycopins_annotate(configuration)

  get_fungal_traits_of_organisms(configuration)

  get_fungi_only_counts(configuration)
  get_wood_saprotrophs_fungi_counts(configuration)

  print("[MycoPins Workflow] Complete!")
}


mycopins_config <- function(batch_name, scata_dataset_name, scata_job_id) {
  reference_data_directory <- "./reference-data/"
  fungal_traits_directory <- paste0(reference_data_directory, "fungaltraits/")
  fungal_traits_file <- paste0(fungal_traits_directory, "fungal_traits_for_genera.csv")

  data_directory <- "./data/"
  batch_directory <- paste0(data_directory, batch_name, "/")

  env_file_name <- paste0(batch_name, "-env.csv")
  env_file <- paste0(batch_directory, env_file_name)

  scata_directory <- paste0(batch_directory, "scata/")
  scata_job_directory <- paste0(scata_directory, scata_job_id, "/")
  scata_counts_file <- paste0(scata_job_directory, "all_tag_by_cluster_counts.txt")
  scata_clusters_file <- paste0(scata_job_directory, "all_clusters_", scata_job_id, ".txt")

  cleaned_counts_file <- paste0(batch_directory, "cleaned_counts.csv")
  cleaned_counts_qc_file <- paste0(batch_directory, "cleaned_counts_qc.txt")
  cleaned_clusters_file <- paste0(batch_directory, "cleaned_clusters.csv")

  blast_directory <- paste0(batch_directory, "blast/")
  consolidated_blast_file <- paste0(batch_directory, "blast_consolidated.csv")
  top_match_file <- paste0(batch_directory, "blast_top_match.csv")

  weather_file <- paste0(batch_directory, "weather.csv")

  complete_clusters_file <- paste0(batch_directory, "complete_clusters.csv")
  complete_counts_file <- paste0(batch_directory, "complete_counts.csv")
  mycopins_environment_file <- paste0(batch_directory, "mycopins_environment.csv")
  mycopins_community_file <- paste0(batch_directory, "mycopins_community.csv")
  mycopins_sample_cluster_file <- paste0(batch_directory, "mycopins_sample_cluster.csv")
  mycopins_organisms_file <- paste0(batch_directory, "mycopins_organisms.csv")
  genus_traits_file <- paste0(batch_directory, "mycopins_genus_traits.csv")

  counts_fungi_file <- paste0(batch_directory, "complete_counts_fungi.csv")
  mycopins_environment_fungi_file <- paste0(batch_directory, "mycopins_environment_fungi.csv")
  mycopins_community_fungi_file <- paste0(batch_directory, "mycopins_community_fungi.csv")

  counts_wood_saprotroph_file <- paste0(
    batch_directory, "complete_counts_wood_saprotrophs.csv")
  mycopins_environment_wood_saprotroph_file <- paste0(
    batch_directory, "mycopins_environment_wood_saprotroph.csv")
  mycopins_community_wood_saprotroph_file <- paste0(
    batch_directory, "mycopins_community_wood_saprotroph.csv")

  return (c(
    data_directory = data_directory,
    batch_directory = batch_directory,
    scata_dataset_name = scata_dataset_name,
    scata_job_id = scata_job_id,
    scata_counts_file = scata_counts_file,
    scata_clusters_file = scata_clusters_file,
    env_file = env_file,
    cleaned_counts_file = cleaned_counts_file,
    cleaned_counts_qc_file = cleaned_counts_qc_file,
    cleaned_clusters_file = cleaned_clusters_file,
    blast_directory = blast_directory,
    consolidated_blast_file = consolidated_blast_file,
    top_match_file = top_match_file,
    weather_file = weather_file,
    complete_clusters_file = complete_clusters_file,
    complete_counts_file = complete_counts_file,
    mycopins_environment_file = mycopins_environment_file,
    mycopins_community_file = mycopins_community_file,
    mycopins_sample_cluster_file = mycopins_sample_cluster_file,
    mycopins_organisms_file = mycopins_organisms_file,
    genus_traits_file = genus_traits_file,
    fungal_traits_file = fungal_traits_file,
    counts_fungi_file = counts_fungi_file,
    mycopins_environment_fungi_file = mycopins_environment_fungi_file,
    mycopins_community_fungi_file = mycopins_community_fungi_file,
    counts_wood_saprotroph_file = counts_wood_saprotroph_file,
    mycopins_environment_wood_saprotroph_file = mycopins_environment_wood_saprotroph_file,
    mycopins_community_wood_saprotroph_file = mycopins_community_wood_saprotroph_file
  ))
}


mycopins_preprocess <- function(configuration) {
  scata_dataset_name = configuration["scata_dataset_name"]
  env_file = configuration["env_file"]
  scata_counts_file = configuration["scata_counts_file"]
  scata_clusters_file = configuration["scata_clusters_file"]
  cleaned_counts_file = configuration["cleaned_counts_file"]
  mycopins_environment_file = configuration["mycopins_environment_file"]
  cleaned_counts_qc_file = configuration["cleaned_counts_qc_file"]
  cleaned_clusters_file = configuration["cleaned_clusters_file"]

  # Clean counts
  print("[Pre-process] Clean counts dataset")
  cleaned_counts_df <- clean_counts(scata_counts_file,
                                    env_file,
                                    scata_dataset_name,
                                    cleaned_counts_qc_file)
  write.csv(cleaned_counts_df, file = cleaned_counts_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  # Clean clusters. clean_clusters writes the data frame to file.
  print("[Pre-process] Clean clusters dataset")
  cleaned_clusters_df <- clean_clusters(scata_clusters_file, cleaned_clusters_file)

  print(paste("[Pre-process] Complete!",
              "Please check 'cleaned_counts_qc.txt' for any ACTION REQUIRED.",
              "DO NOT PROCEED UNTIL ACTION REQUIRED ITEMS ARE RESOLVED."))
}


mycopins_search <- function(configuration) {
  blast_directory = configuration["blast_directory"]
  cleaned_clusters_file = configuration["cleaned_clusters_file"]
  cleaned_clusters_df <- read_csv(cleaned_clusters_file,
                                  show_col_types = FALSE,
                                  locale = locale(encoding = "UTF-8"))

  # Blast cluster sequences
  print("[Search] Create FASTA files")
  batch_cluster_sequences(cleaned_clusters_df, blast_directory)

  print(paste("[Search] FASTA files created.",
    "Upload the FASTA files to BLAST.",
    "Then, download search results in JSON format."))
}


mycopins_identify <- function(configuration) {
  blast_directory = configuration["blast_directory"]
  scata_job_id = configuration["scata_job_id"]
  consolidated_blast_file = configuration["consolidated_blast_file"]
  top_match_file = configuration["top_match_file"]

  # Consolidate BLAST results
  print("[Identify] Consolidate BLAST results")
  consolidated_blast_df <- consolidate_blast_results(blast_directory, scata_job_id)
  write.csv(consolidated_blast_df, file = consolidated_blast_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  # Top BLAST results by cluster
  print("[Identify] Determine top match per cluster")
  top_match_df <- top_blast_results_by_cluster(consolidated_blast_df)
  write.csv(top_match_df, file = top_match_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  print("[Identify] Complete!")
}


mycopins_annotate <- function(configuration) {
  top_match_file = configuration["top_match_file"]
  top_match_df <- read_csv(top_match_file,
                           show_col_types = FALSE,
                           locale = locale(encoding = "UTF-8"))

  cleaned_clusters_file = configuration["cleaned_clusters_file"]
  cleaned_clusters_df <- read_csv(cleaned_clusters_file,
                                  show_col_types = FALSE,
                                  locale = locale(encoding = "UTF-8"))

  cleaned_counts_file = configuration["cleaned_counts_file"]
  cleaned_counts_df <- read_csv(cleaned_counts_file,
                                show_col_types = FALSE,
                                locale = locale(encoding = "UTF-8"))

  weather_file = configuration["weather_file"]
  weather_df <- read_csv(weather_file,
                         show_col_types = FALSE,
                         locale = locale(encoding = "UTF-8"))
  weather_df$date <- as.Date(weather_df$date, "%m/%d/%Y")

  complete_clusters_file = configuration["complete_clusters_file"]
  complete_counts_file = configuration["complete_counts_file"]
  mycopins_environment_file = configuration["mycopins_environment_file"]
  mycopins_community_file = configuration["mycopins_community_file"]
  mycopins_sample_cluster_file = configuration["mycopins_sample_cluster_file"]

  # Apply top match to clusters
  print("[Annotate] Apply top match species to clusters")
  complete_clusters_df <- apply_top_match_to_clusters(top_match_df, cleaned_clusters_df)
  write.csv(complete_clusters_df, file = complete_clusters_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  # Apply match to counts
  print("[Annotate] Apply top match species to counts")
  complete_counts_df <- apply_match_to_counts(cleaned_counts_df, complete_clusters_df,
                                              weather_df, get_weather_data)
  write.csv(complete_counts_df, file = complete_counts_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  sep_index <- which(names(complete_counts_df) == "_Splitter_")
  environment_df <- complete_counts_df[, 1:(sep_index-1)]
  write.csv(environment_df, file = mycopins_environment_file,
            row.names = FALSE, fileEncoding = "UTF-8")
  community_df <- complete_counts_df[, (sep_index+1):ncol(complete_counts_df)]
  write.csv(community_df, file = mycopins_community_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  # Create organism dataset
  print("[Annotate] Create organism dataset")
  organisms <- unique(sort(complete_clusters_df$gbif_accepted_name))

  fungal_traits_file <- configuration["fungal_traits_file"]
  fungal_traits_df <- read_csv(fungal_traits_file,
                               show_col_types = FALSE,
                               locale = locale(encoding = "UTF-8"))

  mycopins_organisms_file <- configuration["mycopins_organisms_file"]
  mycopins_organisms_df <- create_organism_dataset(organisms, complete_clusters_df, fungal_traits_df)
  write.csv(mycopins_organisms_df, file = mycopins_organisms_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  # Create sample-cluster dataset
  print("[Annotate] Create sample-cluster dataset")
  sample_cluster_df <- create_sample_cluster_dataset(cleaned_counts_df, complete_clusters_df, fungal_traits_df)
  write.csv(sample_cluster_df, file = mycopins_sample_cluster_file,
            row.names = FALSE, na = "", fileEncoding = "UTF-8")

  print("[Annotate] Complete!")
}


get_fungal_traits_of_organisms <- function(configuration) {
  print("[Fungal Traits] Retrieve Fungal Traits of cluster organisms")
  mycopins_organisms_file <- configuration["mycopins_organisms_file"]
  mycopins_organisms_df <- read_csv(mycopins_organisms_file,
                                    show_col_types = FALSE,
                                    locale = locale(encoding = "UTF-8"))

  fungal_traits_file <- configuration["fungal_traits_file"]
  fungal_traits_df <- read_csv(fungal_traits_file,
                               show_col_types = FALSE,
                               locale = locale(encoding = "UTF-8"))

  genus_traits_file <- configuration["genus_traits_file"]
  genus_traits_df <- get_fungal_traits(fungal_traits_df, mycopins_organisms_df)
  write.csv(genus_traits_df, file = genus_traits_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  print("[Fungal Traits] Complete")
}


get_fungi_only_counts <- function(configuration) {
  print("[Fungi Only] Filter for fungi only community")
  complete_counts_file = configuration["complete_counts_file"]
  complete_counts_df <- read_csv(complete_counts_file,
                                 show_col_types = FALSE,
                                 locale = locale(encoding = "UTF-8"))

  print("[Fungi Only] Loading organisms dataset")
  mycopins_organisms_file <- configuration["mycopins_organisms_file"]
  mycopins_organisms_df <- read_csv(mycopins_organisms_file,
                                    show_col_types = FALSE,
                                    locale = locale(encoding = "UTF-8"))

  counts_fungi_file <- configuration["counts_fungi_file"]
  counts_fungi_df <- filter_for_fungi(complete_counts_df, mycopins_organisms_df)
  write.csv(counts_fungi_df, file = counts_fungi_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  sep_index <- which(names(counts_fungi_df) == "_Splitter_")
  environment_fungi_df <- counts_fungi_df[, 1:(sep_index-1)]
  mycopins_environment_fungi_file = configuration["mycopins_environment_fungi_file"]
  write.csv(environment_fungi_df, file = mycopins_environment_fungi_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  community_fungi_df <- counts_fungi_df[, (sep_index+1):ncol(counts_fungi_df)]
  mycopins_community_fungi_file = configuration["mycopins_community_fungi_file"]
  write.csv(community_fungi_df, file = mycopins_community_fungi_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  print(paste("[Fungi Only] Complete! Identified",
              ncol(counts_fungi_df), "fungi."))
}


get_wood_saprotrophs_fungi_counts <- function(configuration) {
  print("[Wood Saprotroph Only] Filter for wood saprotroph fungi only community")
  fungal_traits_file <- configuration["fungal_traits_file"]
  fungal_traits_df <- read_csv(fungal_traits_file,
                               show_col_types = FALSE,
                               locale = locale(encoding = "UTF-8"))

  print("[Wood Saprotroph Only] Retrieve wood saprotroph traits")
  wood_saprotrophs_df <- get_wood_saprotrophs_traits(fungal_traits_df)

  print("[Wood Saprotroph Only] Loading complete counts")
  complete_counts_file = configuration["complete_counts_file"]
  complete_counts_df <- read_csv(complete_counts_file,
                                 show_col_types = FALSE,
                                 locale = locale(encoding = "UTF-8"))

  print("[Wood Saprotroph Only] Loading organism dataset")
  mycopins_organisms_file <- configuration["mycopins_organisms_file"]
  mycopins_organisms_df <- read_csv(mycopins_organisms_file,
                                    show_col_types = FALSE,
                                    locale = locale(encoding = "UTF-8"))

  counts_wood_saprotroph_file <- configuration["counts_wood_saprotroph_file"]
  counts_wood_saprotroph_df <- filter_for_wood_saprotroph(
    complete_counts_df, mycopins_organisms_df, wood_saprotrophs_df)

  write.csv(counts_wood_saprotroph_df,file = counts_wood_saprotroph_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  sep_index <- which(names(counts_wood_saprotroph_df) == "_Splitter_")
  environment_wood_saprotroph_df <- counts_wood_saprotroph_df[, 1:(sep_index-1)]
  mycopins_environment_wood_saprotroph_file = configuration["mycopins_environment_wood_saprotroph_file"]
  write.csv(environment_wood_saprotroph_df, file = mycopins_environment_wood_saprotroph_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  community_wood_saprotroph_df <- counts_wood_saprotroph_df[,
                                    (sep_index+1):ncol(counts_wood_saprotroph_df)]
  mycopins_community_wood_saprotroph_file = configuration["mycopins_community_wood_saprotroph_file"]
  write.csv(community_wood_saprotroph_df, file = mycopins_community_wood_saprotroph_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  print(paste("[Wood Saprotroph Only] Complete! Identified",
              ncol(community_wood_saprotroph_df), "wood saprotroph fungi."))
}
