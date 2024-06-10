if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")

library(dplyr)
library(readr)

lib_directory <- paste0(workflow_directory, "lib/")

# Merge
source(paste0(lib_directory, "merge_counts.R"))

mycopins_merge_config <- function(merge_name,
                                  first_name,
                                  second_name) {
  reference_data_directory <- "./reference-data/"

  data_directory <- "./data/"

  merge_directory <- paste0(data_directory, merge_name, "/")
  merge_environment_file <- paste0(merge_directory, "mycopins_environment.csv")
  merge_community_file <- paste0(merge_directory, "mycopins_community.csv")
  merge_environment_fungi_file <- paste0(merge_directory, "mycopins_environment_fungi.csv")
  merge_community_fungi_file <- paste0(merge_directory, "mycopins_community_fungi.csv")
  merge_gbif_taxon_file <- paste0(merge_directory, "mycopins_gbif_taxon.csv")

  first_directory <- paste0(data_directory, first_name, "/")
  first_environment_file <- paste0(first_directory, "mycopins_environment.csv")
  first_community_file <- paste0(first_directory, "mycopins_community.csv")
  first_environment_fungi_file <- paste0(first_directory, "mycopins_environment_fungi.csv")
  first_community_fungi_file <- paste0(first_directory, "mycopins_community_fungi.csv")
  first_gbif_taxon_file <- paste0(first_directory, "mycopins_gbif_taxon.csv")

  second_directory <- paste0(data_directory, second_name, "/")
  second_environment_file <- paste0(second_directory, "mycopins_environment.csv")
  second_community_file <- paste0(second_directory, "mycopins_community.csv")
  second_environment_fungi_file <- paste0(second_directory, "mycopins_environment_fungi.csv")
  second_community_fungi_file <- paste0(second_directory, "mycopins_community_fungi.csv")
  second_gbif_taxon_file <- paste0(second_directory, "mycopins_gbif_taxon.csv")

  return(c(
    merge_directory = merge_directory,
    merge_environment_file = merge_environment_file,
    merge_community_file = merge_community_file,
    merge_environment_fungi_file = merge_environment_fungi_file,
    merge_community_fungi_file = merge_community_fungi_file,
    merge_gbif_taxon_file = merge_gbif_taxon_file,
    first_directory = first_directory,
    first_environment_file = first_environment_file,
    first_community_file = first_community_file,
    first_environment_fungi_file = first_environment_fungi_file,
    first_community_fungi_file = first_community_fungi_file,
    first_gbif_taxon_file = first_gbif_taxon_file,
    second_directory = second_directory,
    second_environment_file = second_environment_file,
    second_community_file = second_community_file,
    second_environment_fungi_file = second_environment_fungi_file,
    second_community_fungi_file = second_community_fungi_file,
    second_gbif_taxon_file = second_gbif_taxon_file
  ))
}

mycopins_merge <- function(configuration) {
  print("[Merge] Merge all counts")
  mycopins_merge_counts(configuration["first_environment_file"],
                        configuration["first_community_file"],
                        configuration["second_environment_file"],
                        configuration["second_community_file"],
                        configuration["merge_environment_file"],
                        configuration["merge_community_file"])
  print("[Merge] Merge fungi only counts")
  mycopins_merge_counts(configuration["first_environment_fungi_file"],
                        configuration["first_community_fungi_file"],
                        configuration["second_environment_fungi_file"],
                        configuration["second_community_fungi_file"],
                        configuration["merge_environment_fungi_file"],
                        configuration["merge_community_fungi_file"])

  print("[Merge] GBIF Taxonomy")
  mycopins_merge_gbif_taxon(configuration["first_gbif_taxon_file"],
                            configuration["second_gbif_taxon_file"],
                            configuration["merge_gbif_taxon_file"])
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


mycopins_merge_gbif_taxon <- function(first_gbif_taxon_file,
                                      second_gbif_taxon_file,
                                      merge_gbif_taxon_file) {
  # Merge environments
  first.gbif_taxon <- read_csv(first_gbif_taxon_file, show_col_types = FALSE)
  second.gbif_taxon <- read_csv(second_gbif_taxon_file, show_col_types = FALSE)
  merged.gbif_taxon <- rbind(first.gbif_taxon, second.gbif_taxon)
  merged.gbif_taxon <- cbind(rowId = rownames(merged.gbif_taxon), merged.gbif_taxon)

  organisms <- unique(subset(merged.gbif_taxon, taxonKey != "KINGDOM")$organism)
  organisms <- setdiff(organisms, c("alga",
                                    "fungi", "fungus",
                                    "Fungi", "Fungus",
                                    "mock", "uncultured organism"))

  rows_to_keep <- c()
  for (organism_name in organisms) {
    search_matches <- merged.gbif_taxon %>%
      filter(organism == organism_name & searchBitScore >= 200) %>%
      arrange(desc(searchPercentIdentity), searchEValue)

    best_search_match <- NULL
    if (nrow(search_matches) > 0) {
      rows_to_keep <- c(rows_to_keep, search_matches[1,]$rowId)
    }
  }

  # Keep only best search match rows for each organism
  merged.gbif_taxon <- subset(merged.gbif_taxon, rowId %in% rows_to_keep)
  # Remove rowId column
  merged.gbif_taxon <- subset(merged.gbif_taxon, select = -rowId)
  # Sort based on organism
  merged.gbif_taxon <- merged.gbif_taxon[order(merged.gbif_taxon$organism), ]

  write.csv(merged.gbif_taxon, file = merge_gbif_taxon_file, row.names = FALSE)
}
