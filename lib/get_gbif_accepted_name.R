if (!require("rgbif")) install.packages("rgbif")

library("rgbif")

get_gbif_accepted_name <- function(normalized_name) {
  accepted_name <- normalized_name

  response <- name_backbone_verbose(name=normalized_name, kingdom="fungi")
  taxon <- response$data

  accepted_name_usage <- ifelse("species" %in% names(taxon), taxon$species, "")
  original_name_usage <- ifelse("canonicalName" %in% names(taxon),
                                taxon$canonicalName, "")

  if (taxon$matchType != "NONE") {
    if (!is.na(accepted_name_usage) & nchar(accepted_name_usage) > 0) {
      accepted_name <- accepted_name_usage
    } else if (!is.na(original_name_usage) & nchar(original_name_usage) > 0) {
      accepted_name <- original_name_usage
    }
  }

  return(accepted_name)
}