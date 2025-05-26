lib_directory <- paste0(workflow_directory, "lib/")

source(paste0(lib_directory, "generate_gbif_event_core.R"))
source(paste0(lib_directory, "generate_gbif_event_occurrence.R"))
source(paste0(lib_directory, "generate_gbif_dna_derived.R"))
source(paste0(lib_directory, "generate_gbif_occurrence_core.R"))

#
# Helper functions
#
get_transects <- function() {
  return(
    c("A", "B", "C")
  )
}

get_habitat <- function(transect) {
  habitats <- c(
    "A" = "buried in soil of a conifer forest area with reindeer access.",
    "B" = "buried in soil of a conifer forest area protected from grazing by reindeers.",
    "C" = "buried in soil of a mixed broadleaf forest accessed by random visitors."
  )

  return(habitats[transect])
}

get_latitude <- function(transect) {
  latitudes <- c(
    "A" = 66.367,
    "B" = 66.367, # A and B have the same coordinates
    "C" = 66.376
  )

  return(latitudes[transect])
}

get_longitude <- function(transect) {
  longitudes <- c(
    "A" = 29.533,
    "B" = 29.533, # A and B have the same coordinates
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
  return("UNITE Fungi 9.0 (2023-07-18) (https://doi.org/10.15156/BIO/2938068)")
}

#
# Main functions: mycopins_ipt_gbif_config(), mycopins_ipt_gbif_generate(configuration)
#
mycopins_ipt_gbif_config <- function() {
  ipt_gbif_data_directory <- "./data/All"
  ipt_gbif_output_directory <- "./data/IPT-GBIF"

  mycopins_environment_file <- paste0(ipt_gbif_data_directory, "/mycopins_environment.csv")
  mycopins_community_file <- paste0(ipt_gbif_data_directory, "/mycopins_community.csv")
  mycopins_sample_cluster_file <- paste0(ipt_gbif_data_directory, "/mycopins_sample_cluster.csv")
  mycopins_organisms_file <- paste0(ipt_gbif_data_directory, "/mycopins_organisms.csv")

  ipt_gbif_event_core_file <- paste0(ipt_gbif_output_directory, "/event.txt")
  ipt_gbif_event_occurrence_all_fungi_file <- paste0(ipt_gbif_output_directory, "/occurrence-all-fungi.txt")
  ipt_gbif_event_occurrence_file <- paste0(ipt_gbif_output_directory, "/occurrence.txt")
  ipt_gbif_dna_derived_file <- paste0(ipt_gbif_output_directory, "/dna-derived.txt")
  ipt_gbif_occurrence_core_file <- paste0(ipt_gbif_output_directory, "/occurrence-core.txt")

  #
  # Constant definitions
  #
  location <- "Oulanka"

  # Event
  samplingProtocol <- "https://doi.org/10.3897/mycokeys.96.101033"
  locationID <- "https://www.geonames.org/12226273"
  countryCode <- "FI"
  country <- "Finland"
  stateProvince <- "North Ostrobothnia"
  municipality <- "Kuusamo"

  # Occurrence
  basisOfRecord <- "MaterialSample"
  organismQuantityType <- "DNA sequence reads"

  # DNA-derived
  # https://rs.gbif.org/extension/gbif/1.0/dna_derived_data_2024-04-17.xml#DNA_sequence
  libLayout <- "paired"
  targetGene <- "ITS"
  targetSubfragment <- "ITS2"
  seqMethod <- "Illumina MiSeq"
  otuDB <- "Clustering based on the UNITE Fungi 9.0 (2023-07-18) (https://doi.org/10.15156/BIO/2938068) using USearch with 90% identity parameter (SCATA)"
  pcrPrimerNameFormward <- "fITS7"
  pcrPrimerNameReverse <- "ITS4"
  pcrPrimerReference <- "https://doi.org/10.1007/978-1-4939-3369-3_4"

  return(c(
    ipt_gbif_output_directory = ipt_gbif_output_directory,
    ipt_gbif_event_core_file = ipt_gbif_event_core_file,
    ipt_gbif_event_occurrence_all_fungi_file = ipt_gbif_event_occurrence_all_fungi_file,
    ipt_gbif_event_occurrence_file = ipt_gbif_event_occurrence_file,
    ipt_gbif_dna_derived_file = ipt_gbif_dna_derived_file,
    ipt_gbif_occurrence_core_file = ipt_gbif_occurrence_core_file,
    mycopins_environment_file = mycopins_environment_file,
    mycopins_community_file = mycopins_community_file,
    mycopins_sample_cluster_file = mycopins_sample_cluster_file,
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
    # DNA-derived data
    libLayout = libLayout,
    targetGene = targetGene,
    targetSubfragment = targetSubfragment,
    seqMethod = seqMethod,
    otuDB = otuDB,
    pcrPrimerNameFormward = pcrPrimerNameFormward,
    pcrPrimerNameReverse = pcrPrimerNameReverse,
    pcrPrimerReference = pcrPrimerReference
  ))
}


mycopins_ipt_gbif_generate <- function(configuration) {
  #
  # event core
  #
  # gbif_event_core <- generate_gbif_event_core(configuration)
  # write.csv(gbif_event_core, configuration["ipt_gbif_event_core_file"],
  #     row.names = FALSE, na = "", fileEncoding = "UTF-8")


  #
  # event occurrence (extension)
  #
  gbif_event_occurrence <- generate_gbif_event_occurrence(configuration)
  write.csv(gbif_event_occurrence, configuration["ipt_gbif_event_occurrence_all_fungi_file"],
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

  gbif_event_occurrence_no_kingdom <- generate_gbif_event_occurrence_no_kingdom(configuration)
  write.csv(gbif_event_occurrence_no_kingdom, configuration["ipt_gbif_event_occurrence_file"],
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

  #
  # dna-derived
  #
  # gbif_dna_derived <- generate_dna_derived(configuration)
  # write.csv(gbif_dna_derived, configuration["ipt_gbif_dna_derived_file"],
  #         row.names = FALSE, na = "", fileEncoding = "UTF-8")


  #
  # occurrence core
  #
  # gbif_occurrence_core <- generate_gbif_occurrence_core(configuration)
  # write.csv(gbif_occurrence_core, configuration["ipt_gbif_occurrence_core_file"],
  #         row.names = FALSE, na = "", fileEncoding = "UTF-8")
}