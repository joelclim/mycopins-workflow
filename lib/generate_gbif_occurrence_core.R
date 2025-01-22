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

  occurrence_core_df <- occurrence_core_df %>%
    mutate(
      dynamicProperties = mapply(function(x, y) {
        json1 <- fromJSON(x)
        json2 <- fromJSON(y)
        json2$woodType <- NULL
        json2$woodTexture <- NULL

        # Merge the two JSON objects
        merged_json <- toJSON(c(json1, json2), auto_unbox = TRUE)
        return(merged_json)
      }, dynamicProperties.x, dynamicProperties.y)
    ) %>%
    select(-dynamicProperties.x, -dynamicProperties.y)

  return(occurrence_core_df)
}
