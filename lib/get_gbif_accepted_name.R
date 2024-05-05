library("rgbif")

get_gbif_accepted_name <- function(normalized_name) {
  accepted_name <- normalized_name
  
  response <- name_backbone_verbose(name=normalized_name, kingdom="fungi")
  taxon <- response$data

  acceptedNameUsage <- ifelse("species" %in% names(taxon), taxon$species, "")
  originalNameUsage <- ifelse("canonicalName" %in% names(taxon), taxon$canonicalName, "")
  
  if (taxon$matchType != "NONE") { 
    if (!is.na(acceptedNameUsage) & nchar(acceptedNameUsage) > 0) {
      accepted_name <- acceptedNameUsage
    } else if (!is.na(originalNameUsage) & nchar(originalNameUsage) > 0) {
      accepted_name <- originalNameUsage
    }
  }
  
  print(paste0("GBIF accepted name for '", normalized_name, 
               "' is '", accepted_name, "'"))
  
  return(accepted_name) 
}