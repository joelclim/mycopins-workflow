if (!require("readr")) install.packages("readr")

library(readr)

lib_directory <- paste0(workflow_directory, "lib/")

# Merge
source(paste0(lib_directory, "merge_counts.R"))

# Taxonomy and Traits
source(paste0(lib_directory, "get_gbif_taxon.R"))
source(paste0(lib_directory, "get_fungal_traits.R"))

# Filter for Fungi
source(paste0(lib_directory, "filter_for_fungi.R"))
source(paste0(lib_directory, "filter_for_wood_saprotroph.R"))

mycopins_merge_config <- function(merge_name,
                                  first_name,
                                  second_name) {
  data_directory <- "./data/"

  merge_directory <- paste0(data_directory, merge_name, "/")
  merge_environment_file <- paste0(merge_directory, "mycopins_environment.csv")
  merge_community_file <- paste0(merge_directory, "mycopins_community.csv")
  merge_environment_fungi_file <- paste0(merge_directory, "mycopins_environment_fungi.csv")
  merge_community_fungi_file <- paste0(merge_directory, "mycopins_community_fungi.csv")
  merge_environment_wood_saprotroph_file <- paste0(merge_directory,
                                                 "mycopins_environment_wood_saprotroph.csv")
  merge_community_wood_saprotroph_file <- paste0(merge_directory,
                                                 "mycopins_community_wood_saprotroph.csv")

  first_directory <- paste0(data_directory, first_name, "/")
  first_environment_file <- paste0(first_directory, "mycopins_environment.csv")
  first_community_file <- paste0(first_directory, "mycopins_community.csv")
  first_environment_fungi_file <- paste0(first_directory, "mycopins_environment_fungi.csv")
  first_community_fungi_file <- paste0(first_directory, "mycopins_community_fungi.csv")
  first_environment_wood_saprotroph_file <- paste0(first_directory,
                                                 "mycopins_environment_wood_saprotroph.csv")
  first_community_wood_saprotroph_file <- paste0(first_directory,
                                                 "mycopins_community_wood_saprotroph.csv")

  second_directory <- paste0(data_directory, second_name, "/")
  second_environment_file <- paste0(second_directory, "mycopins_environment.csv")
  second_community_file <- paste0(second_directory, "mycopins_community.csv")
  second_environment_fungi_file <- paste0(second_directory, "mycopins_environment_fungi.csv")
  second_community_fungi_file <- paste0(second_directory, "mycopins_community_fungi.csv")
  second_environment_wood_saprotroph_file <- paste0(second_directory,
                                                 "mycopins_environment_wood_saprotroph.csv")
  second_community_wood_saprotroph_file <- paste0(second_directory,
                                                 "mycopins_community_wood_saprotroph.csv")

  return(c(
    merge_directory = merge_directory,
    merge_environment_file = merge_environment_file,
    merge_community_file = merge_community_file,
    merge_environment_fungi_file = merge_environment_fungi_file,
    merge_community_fungi_file = merge_community_fungi_file,
    merge_environment_wood_saprotroph_file = merge_environment_wood_saprotroph_file,
    merge_community_wood_saprotroph_file = merge_community_wood_saprotroph_file,
    first_directory = first_directory,
    first_environment_file = first_environment_file,
    first_community_file = first_community_file,
    first_environment_fungi_file = first_environment_fungi_file,
    first_community_fungi_file = first_community_fungi_file,
    first_environment_wood_saprotroph_file = first_environment_wood_saprotroph_file,
    first_community_wood_saprotroph_file = first_community_wood_saprotroph_file,
    second_directory = second_directory,
    second_environment_file = second_environment_file,
    second_community_file = second_community_file,
    second_environment_fungi_file = second_environment_fungi_file,
    second_community_fungi_file = second_community_fungi_file,
    second_environment_wood_saprotroph_file = second_environment_wood_saprotroph_file,
    second_community_wood_saprotroph_file = second_community_wood_saprotroph_file
  ))
}

mycopins_merge <- function(configuration) {
  print("[Merge] All counts")
  mycopins_merge_counts(configuration["first_environment_file"],
                        configuration["first_community_file"],
                        configuration["second_environment_file"],
                        configuration["second_community_file"],
                        configuration["merge_environment_file"],
                        configuration["merge_community_file"])
  print("[Merge] Fungi only counts")
  mycopins_merge_counts(configuration["first_environment_fungi_file"],
                        configuration["first_community_fungi_file"],
                        configuration["second_environment_fungi_file"],
                        configuration["second_community_fungi_file"],
                        configuration["merge_environment_fungi_file"],
                        configuration["merge_community_fungi_file"])
  print("[Merge] Wood Saprotrophs only counts")
  mycopins_merge_counts(configuration["first_environment_wood_saprotroph_file"],
                        configuration["first_community_wood_saprotroph_file"],
                        configuration["second_environment_wood_saprotroph_file"],
                        configuration["second_community_wood_saprotroph_file"],
                        configuration["merge_environment_wood_saprotroph_file"],
                        configuration["merge_community_wood_saprotroph_file"])
}

mycopins_merge_counts <- function(first_environment_file,
                                  first_community_file,
                                  second_environment_file,
                                  second_community_file,
                                  merge_environment_file,
                                  merge_community_file) {
  # Merge environments
  first.env <- read_csv(first_environment_file, show_col_types = FALSE)
  second.env <- read_csv(second_environment_file, show_col_types = FALSE)
  merged.env <- rbind(first.env, second.env)
  write.csv(merged.env, file = merge_environment_file, row.names = FALSE)

  # Merge communities
  first <- read_csv(first_community_file, show_col_types = FALSE)
  second <- read_csv(second_community_file, show_col_types = FALSE)
  merged <- merge_communities(merged.env, first, second)
  write.csv(merged, file = merge_community_file, row.names = FALSE)
}
