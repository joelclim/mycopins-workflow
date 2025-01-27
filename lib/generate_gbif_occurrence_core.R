if (!require("dplyr")) install.packages("dplyr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("readr")) install.packages("readr")

library(dplyr)
library(readr)
library(jsonlite)

#
# Main function: generate_dna_derived
#
generate_gbif_occurrence_core <- function(configuration) {
  event_core_file <- configuration["ipt_gbif_event_core_file"]
  event_core_df <- read_csv(event_core_file, show_col_types = FALSE,
                                locale = locale(encoding = "UTF-8"))

  event_occurrence_file <- configuration["ipt_gbif_event_occurrence_file"]
  event_occurrence_df <- read_csv(event_occurrence_file, show_col_types = FALSE,
                            locale = locale(encoding = "UTF-8"))

  occurrence_core_df <- event_occurrence_df %>%
                        inner_join(event_core_df, by = "eventID")

  merge_dynamic_properties <- function(json_str1, json_str2) {
    list1 <- fromJSON(json_str1)
    list2 <- fromJSON(json_str2)

    # Combine the lists
    merged_list <- list1

    # Iterate through the keys of list2
    for (key in names(list2)) {
      if (!key %in% names(list1)) {
        merged_list[[key]] <- list2[[key]]
      }
    }

    # Convert the merged list back to JSON
    toJSON(merged_list, auto_unbox = TRUE)
  }

  occurrence_core_df$dynamicProperties <- mapply(merge_dynamic_properties,
                          occurrence_core_df$dynamicProperties.x,
                          occurrence_core_df$dynamicProperties.y)

  occurrence_core_df <- occurrence_core_df[,
    !(names(occurrence_core_df) %in% c("dynamicProperties.x", "dynamicProperties.y"))]

  return(occurrence_core_df)
}
