lib_directory <- paste0(workflow_directory, "lib/")

source(paste0(lib_directory, "generate_gbif_event.R"))
source(paste0(lib_directory, "generate_gbif_occurrence.R"))

#
# Helper functions
#
get_transects <- function() {
  return(
    c("A", "C")
  )
}

get_habitat <- function(transect) {
  habitats <- c(
    "A" = "buried in soil of a boreal forest area protected from grazing by reindeers.",
    "C" = "buried in soil of a mixed broadleaf forest accessed by random visitors."
  )

  return(habitats[transect])
}

get_latitude <- function(transect) {
  latitudes <- c(
    "A" = 66.367,
    "C" = 66.376
  )

  return(latitudes[transect])
}

get_longitude <- function(transect) {
  longitudes <- c(
    "A" = 29.533,
    "C" = 29.313
  )

  return(longitudes[transect])
}

get_wood_type <- function(site) {
  if (site %in% c("A", "B")) {
    return("Pine")
  }
  if (site %in% c("C", "D")) {
    return("Birch")
  }
  # site %in% c("E", "F")
  return ("Spruce")
}

get_wood_texture <- function(site) {
  if (site %in% c("A", "B", "E", "F")) {
    return("Softwood")
  }
  # site %in% c("C", "D")
  return ("Hardwood")
}

get_match_source <- function(source) {
  if (source == "BLAST") {
    return("BLASTn circa 2024-05-01")
  }
  return("UNITE Fungi 9.0 (2023-07-18)")
}

#
# Main functions: mycopins_ipt_gbif_config(), mycopins_ipt_gbif_generate(configuration)
#
mycopins_ipt_gbif_config <- function() {
  ipt_gbif_data_directory <- "./data/All"
  ipt_gbif_output_directory <- "./data/IPT-GBIF"

  mycopins_environment_file <- paste0(ipt_gbif_data_directory, "/mycopins_environment_fungi.csv")
  mycopins_community_file <- paste0(ipt_gbif_data_directory, "/mycopins_community_fungi.csv")
  mycopins_organisms_file <- paste0(ipt_gbif_data_directory, "/mycopins_organisms_fungi.csv")

  ipt_gbif_event_file <- paste0(ipt_gbif_output_directory, "/event.txt")
  ipt_gbif_occurrence_file <- paste0(ipt_gbif_output_directory, "/occurrence.txt")

  #
  # Constant definitions
  #
  location <- "Oulanka"

  # Event
  samplingProtocol <- "Shumskaya, M., Lorusso, N., Patel, U., Leigh, M., Somervuo, P., & Schigel, D. (2023). MycoPins: a metabarcoding-based method to monitor fungal colonization of fine woody debris. MycoKeys, 96, 77–95. https://doi.org/10.3897/mycokeys.96.101033"
  locationID <- "https://www.geonames.org/12226273"
  countryCode <- "FI"
  country <- "Finland"
  stateProvince <- "North Ostrobothnia"
  municipality <- "Kuusamo"

  # Occurrence
  basisOfRecord <- "MaterialSample"
  organismQuantityType <- "DNA sequence reads"
  sampleSizeUnit <- "DNA sequence reads"

  return(c(
    ipt_gbif_output_directory = ipt_gbif_output_directory,
    ipt_gbif_event_file = ipt_gbif_event_file,
    ipt_gbif_occurrence_file = ipt_gbif_occurrence_file,
    mycopins_environment_file = mycopins_environment_file,
    mycopins_community_file = mycopins_community_file,
    mycopins_organisms_file = mycopins_organisms_file,
    # Shared constants
    location = location,
    # Event
    samplingProtocol  = samplingProtocol,
    locationID = locationID,
    countryCode = countryCode,
    country = country,
    stateProvince = stateProvince,
    municipality = municipality,
    # Occurrence
    basisOfRecord = basisOfRecord,
    organismQuantityType = organismQuantityType,
    sampleSizeUnit = sampleSizeUnit
  ))
}


mycopins_ipt_gbif_generate <- function(configuration) {
  #
  # event
  #
  gbif_event <- generate_gbif_event(configuration)
  write.csv(gbif_event, configuration["ipt_gbif_event_file"],
      row.names = FALSE, na = "", fileEncoding = "UTF-8")


  #
  # occurrence
  #
  gbif_occurrence <- generate_gbif_occurrence(configuration)
  write.csv(gbif_occurrence, configuration["ipt_gbif_occurrence_file"],
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
}