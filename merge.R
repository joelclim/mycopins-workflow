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
  merge_organisms_file <- paste0(merge_directory, "mycopins_organisms.csv")
  merge_organisms_fungi_file <- paste0(merge_directory, "mycopins_organisms_fungi.csv")

  first_directory <- paste0(data_directory, first_name, "/")
  first_environment_file <- paste0(first_directory, "mycopins_environment.csv")
  first_community_file <- paste0(first_directory, "mycopins_community.csv")
  first_environment_fungi_file <- paste0(first_directory, "mycopins_environment_fungi.csv")
  first_community_fungi_file <- paste0(first_directory, "mycopins_community_fungi.csv")
  first_organisms_file <- paste0(first_directory, "mycopins_organisms.csv")

  second_directory <- paste0(data_directory, second_name, "/")
  second_environment_file <- paste0(second_directory, "mycopins_environment.csv")
  second_community_file <- paste0(second_directory, "mycopins_community.csv")
  second_environment_fungi_file <- paste0(second_directory, "mycopins_environment_fungi.csv")
  second_community_fungi_file <- paste0(second_directory, "mycopins_community_fungi.csv")
  second_organisms_file <- paste0(second_directory, "mycopins_organisms.csv")

  return(c(
    merge_directory = merge_directory,
    merge_environment_file = merge_environment_file,
    merge_community_file = merge_community_file,
    merge_environment_fungi_file = merge_environment_fungi_file,
    merge_community_fungi_file = merge_community_fungi_file,
    merge_organisms_file = merge_organisms_file,
    merge_organisms_fungi_file = merge_organisms_fungi_file,
    first_directory = first_directory,
    first_environment_file = first_environment_file,
    first_community_file = first_community_file,
    first_environment_fungi_file = first_environment_fungi_file,
    first_community_fungi_file = first_community_fungi_file,
    first_organisms_file = first_organisms_file,
    second_directory = second_directory,
    second_environment_file = second_environment_file,
    second_community_file = second_community_file,
    second_environment_fungi_file = second_environment_fungi_file,
    second_community_fungi_file = second_community_fungi_file,
    second_organisms_file = second_organisms_file
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

  print("[Merge] Organisms datasets")
  mycopins_merge_organisms(configuration["first_organisms_file"],
                           configuration["second_organisms_file"],
                           configuration["merge_organisms_file"])

  mycopins_merge_organisms(configuration["first_organisms_file"],
                           configuration["second_organisms_file"],
                           configuration["merge_organisms_fungi_file"], TRUE)
}


mycopins_merge_counts <- function(first_environment_file,
                                  first_community_file,
                                  second_environment_file,
                                  second_community_file,
                                  merge_environment_file,
                                  merge_community_file) {
  # Merge environments
  first.env <- read_csv(first_environment_file,
                        show_col_types = FALSE,
                        locale = locale(encoding = "UTF-8"))
  second.env <- read_csv(second_environment_file,
                         show_col_types = FALSE,
                         locale = locale(encoding = "UTF-8"))
  merged.env <- rbind(first.env, second.env)
  write.csv(merged.env, file = merge_environment_file,
            row.names = FALSE, fileEncoding = "UTF-8")

  # Merge communities
  first <- read_csv(first_community_file,
                    show_col_types = FALSE,
                    locale = locale(encoding = "UTF-8"))
  second <- read_csv(second_community_file,
                     show_col_types = FALSE,
                     locale = locale(encoding = "UTF-8"))
  merged <- merge_communities(merged.env, first, second)
  write.csv(merged, file = merge_community_file,
            row.names = FALSE, fileEncoding = "UTF-8")
}


mycopins_merge_organisms <- function(first_organisms_file,
                                     second_organisms_file,
                                     merge_organisms_file,
                                     fungi_only = FALSE) {
  # Merge environments
  first.organisms <- read_csv(first_organisms_file,
                              show_col_types = FALSE,
                              locale = locale(encoding = "UTF-8"))
  second.organisms <- read_csv(second_organisms_file,
                               show_col_types = FALSE,
                               locale = locale(encoding = "UTF-8"))
  merged.organisms <- rbind(first.organisms, second.organisms)
  merged.organisms <- cbind(rowId = rownames(merged.organisms), merged.organisms)

  organisms <- unique(merged.organisms$organism)
  if (fungi_only) {
    organisms <- unique(subset(merged.organisms, gbif.taxon_rank != "KINGDOM")$organism)
    organisms <- setdiff(organisms, c("alga",
                                      "fungi", "fungus",
                                      "Fungi", "Fungus",
                                      "mock", "uncultured organism"))
  }

  rows_to_keep <- c()
  for (organism_name in organisms) {
    search_matches <- merged.organisms %>%
      filter(organism == organism_name & match.bit_score >= 200) %>%
      arrange(desc(match.percent_identity), match.evalue)

    best_search_match <- NULL
    if (nrow(search_matches) > 0) {
      rows_to_keep <- c(rows_to_keep, search_matches[1,]$rowId)
    }
  }

  # Keep only best search match rows for each organism
  merged.organisms <- subset(merged.organisms, rowId %in% rows_to_keep)
  # Remove rowId column
  merged.organisms <- subset(merged.organisms, select = -rowId)
  # Sort based on organism
  merged.organisms <- merged.organisms[order(merged.organisms$organism), ]

  write.csv(merged.organisms, file = merge_organisms_file,
            row.names = FALSE, fileEncoding = "UTF-8")
}
