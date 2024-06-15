library("rgbif")

organism <- "Aureobasidium pullulans"
organism <- "Acephala"
#organism <- "Acephala applanata"
gbif_match <- name_lookup(query=organism,
                          status="ACCEPTED",
                          origin="SOURCE",
                          hl=FALSE, 
                          limit=1)
taxon <- gbif_match$data
"key" %in% names(taxon)
taxon$key

organism <- "Aureobasidium pullulans"
organism <- "fungus"
organism <- "Cystobasidium pinicola"
organism <- "Kneiffiella flavipora"
organism <- "Winnie the Pooh"
organism <- "Collophora"
organism <- "Spicatispora fennica"
organism <- "Ceraceomyces tessulatus"
#organism <- "Acephala"
organism <- "Acephala applanata"
result <- name_backbone_verbose(name=organism)
data <- result$data
data$matchType
data$usageKey

# Find the position of the given string
record <- data.frame(list(
  scientificName = data$scientificName,
  acceptedNameUsage = data$species,
  originalNameUsage = data$canonicalName,
  kingdom = data$kingdom,
  phylum = data$phylum,
  class = data$class,
  order = data$order,
  family = data$family,
  genus = data$genus,
  taxonKey = data$rank,
  status = data$status,
  scientificNameAuthorship = trimws(substring(data$scientificName, 
                                              nchar(data$canonicalName)+1,
                                              nchar(data$scientificName))) 
))
record


result <- occ_search(scientificName = organism)
result$data$occurrenceID[1]
result$data$occurrenceStatus[1]


taxonKey <- 2525226
# Get occurrence data for the taxonKey
occurrence <- occ_data(taxonKey = taxonKey, limit = 1)

# Extract image URL from metadata
image_url <- occurrence$data$media$type == "StillImage"
if (any(image_url)) {
  image_url <- occurrence$data$media$url[image_url][1]
} else {
  image_url <- NA
}
occurrence$data$
