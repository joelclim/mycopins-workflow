library(readr)

lib_directory <- paste0(workflow_directory, "lib/")

# Merge
source(paste0(lib_directory, "merge_add_counts.R"))

# Taxonomy and Traits
source(paste0(lib_directory, "get_gbif_taxon.R"))
source(paste0(lib_directory, "get_fungal_traits.R"))

# Filter for Fungi
source(paste0(lib_directory, "filter_for_fungi.R"))
source(paste0(lib_directory, "filter_for_wood_saprotroph.R"))

mycopins_merge_config <- function(merge_name,
                                  first_name,
                                  second_name,
                                  configuration) {
  data_directory <- "./data/"

  merge_directory <- paste0(data_directory, merge_name, "/")
  merge_counts_file <- paste0(merge_directory, "merged_counts.csv")

  first_directory <- paste0(data_directory, first_name, "/")
  first_counts_file <- paste0(first_directory, "cleaned_counts.csv")

  second_directory <- paste0(data_directory, second_name, "/")
  second_counts_file <- paste0(second_directory, "cleaned_counts.csv")

  return(c(configuration,
    merge_directory = merge_directory,
    merge_counts_file = merge_counts_file,
    first_directory = first_directory,
    first_counts_file = first_counts_file,
    second_directory = second_directory,
    second_counts_file = second_counts_file
  ))
}

mycopins_merge <- function(merge_configuration) {
  first_counts_file <- merge_configuration["first_counts_file"]
  first_df <- read_csv(first_counts_file)
  first_df <- cbind(first_df)

  second_counts_file <- merge_configuration["second_counts_file"]
  second_df <- read_csv(second_counts_file)
  second_df <- cbind(second_df)

  merge_counts_file <- merge_configuration["merge_counts_file"]
  merged_df <- merge_add_counts(first_df, second_df)
  write.csv(merged_df, file = merge_counts_file, row.names = FALSE)
}