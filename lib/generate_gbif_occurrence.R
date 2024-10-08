library(dplyr)
library(jsonlite)
library(readr)

#
# Requires get_transects(), get_habitat(transect), get_latitude(transect), get_longitude(transect)
#          get_wood_type(site_letter), get_wood_texture(site_letter), get_match_source(source)
#

create_gbif_occurrence_dynamic_properties <- function(site_letter, organism) {
  dynamic_properties_df <- data.frame(woodType = get_wood_type(site_letter),
                                      woodTexture = get_wood_texture(site_letter))
  match_df <- data.frame(
    source = get_match_source(organism$match.source),
    percentIdentity = organism$match.percent_identity,
    evalue = organism$match.evalue,
    bitScore = organism$match.bit_score,
    accession = organism$match.accession,
    taxId = organism$match.tax_id
  )
  dynamic_properties_df$match <- match_df

  fungal_traits_df <- data.frame(
    primaryLifestyle = organism$fungal_traits.primary_lifestyle,
    secondaryLifestyle = organism$fungal_traits.secondary_lifestyle
  )
  dynamic_properties_df$fungalTraits <- fungal_traits_df

  return(dynamic_properties_df)
}

create_gbif_occurrence_record <- function(eventID, occurrenceID,
                                    basisOfRecord, organismQuantityType, sampleSizeUnit,
                                    site_letter, specie, organism, occurrence_match) {
  get_data_generalizations <- function(organism) {
    if (organism$match.source == "BLAST") {
      return(paste0("Identified via BLASTn circa 2024-05-01 with ",
                  organism$match.bit_score, " bit score, ",
                  organism$match.evalue, " evalue, and ",
                  format(round(organism$match.percent_identity, 2), nsmall = 2), "% identity."))
    }
    return("Identified via UNITE Fungi 9.0 (2023-07-18).")
  }

  get_taxon_id <- function(organism) {
    return(paste0("https://www.gbif.org/species/", organism$gbif.usage_key))
  }

  get_taxon_concept_id <- function(organism) {
    if (organism$match.source == "BLAST") {
      return(paste0("NCBI:tx", organism$match.tax_id))
    }

    # UNITE
    tax_id <- strsplit(organism$match.tax_id, split="\\|")
    sh_id <- tax_id[[1]][1]
    seq_dataset <- tax_id[[1]][2]
    return(paste0("UNITE:", sh_id))
  }


  taxonID <- get_taxon_id(organism)
  taxonConceptID <- get_taxon_concept_id(organism)
  dataGeneralizations <- get_data_generalizations(organism)

  organismQuantity <- 0
  if (nrow(occurrence_match) > 0) {
    organismQuantity <- occurrence_match[[ specie ]]
  }

  occurrenceStatus <- "absent"
  if (organismQuantity > 0) {
    occurrenceStatus <- "present"
  }

  dynamic_properties_df <- create_gbif_occurrence_dynamic_properties(site_letter, organism)

  preparations <- paste("DNA extract from a", get_wood_type(site_letter), "MycoPin.")

  occurrence_record <- list(
    occurrenceID = occurrenceID,
    eventID = eventID,
    basisOfRecord = basisOfRecord,
    dataGeneralizations = dataGeneralizations,
    taxonID = taxonID,
    taxonConceptID = taxonConceptID,
    organismQuantity = organismQuantity,
    organismQuantityType = organismQuantityType,
    occurrenceStatus = occurrenceStatus,
    preparations = preparations,
    sampleSizeUnit = sampleSizeUnit,
    scientificName = organism$gbif.scientific_name,
    acceptedNameUsage = organism$gbif.accepted_name_usage,
    originalNameUsage = organism$gbif.original_name_usage,
    kingdom = organism$gbif.kingdom,
    phylum = organism$gbif.phylum,
    class = organism$gbif.class,
    order = organism$gbif.order,
    family = organism$gbif.family,
    genus = organism$gbif.genus,
    taxonRank = organism$gbif.taxon_rank,
    scientificNameAuthorship = organism$gbif.scientific_name_authorship,
    taxonomicStatus = organism$gbif.status,
    dynamicProperties = gsub("\\[|\\]", "", toJSON(dynamic_properties_df))
  )

  return(occurrence_record)
}

#
# Main function: generate_gbif_occurrence
#
generate_gbif_occurrence <- function(configuration) {
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }

  location <- configuration["location"]
  basisOfRecord <- configuration["basisOfRecord"]
  organismQuantityType <- configuration["organismQuantityType"]
  sampleSizeUnit <- configuration["sampleSizeUnit"]

  transects <- get_transects()

  environment_file <- configuration["mycopins_environment_file"]
  mycopins.environment <- read_csv(environment_file, show_col_types = FALSE,
                                  locale = locale(encoding = "UTF-8"))

  community_file <- configuration["mycopins_community_file"]
  mycopins.community <- read_csv(community_file, show_col_types = FALSE,
                                locale = locale(encoding = "UTF-8"))
  mycopins <- cbind(mycopins.environment, mycopins.community)

  organisms_file <- configuration["mycopins_organisms_file"]
  mycopins.organisms <- read_csv(organisms_file, show_col_types = FALSE,
                                locale = locale(encoding = "UTF-8"))

  gbif_occurrences <- data.frame(
    occurrenceID = character(),
    eventID = character(),
    basisOfRecord = character(),
    dataGeneralizations = character(),
    taxonID = character(),
    taxonConceptID = character(),
    organismQuantity = numeric(),
    organismQuantityType = character(),
    occurrenceStatus = character(),
    preparations = character(),
    sampleSizeUnit = character(),
    scientificName = character(),
    acceptedNameUsage = character(),
    originalNameUsage = character(),
    kingdom = character(),
    phylum = character(),
    class = character(),
    order = character(),
    family = character(),
    genus = character(),
    taxonRank = character(),
    scientificNameAuthorship = character(),
    taxonomicStatus = character()
  )

  species <- unique(sort(names(mycopins.community)))

  execution_time <- system.time({
    transect_count <- length(transects)
    for (i in 1:transect_count) {
      transect <- transects[i]
      collection_dates <- mycopins.environment$Date.Collected[
        mycopins.environment$Transect == transect
      ]
      collection_dates <- as.Date(collection_dates, "%m/%d/%Y")
      collection_dates <- sort(unique(collection_dates))

      collection_date_count <- length(collection_dates)
      for (j in 1:collection_date_count) {
        eventDate <- collection_dates[j]
        date_id <- format(eventDate, "%Y_%b_%d")

        sites <- mycopins.environment %>%
          filter(as.Date(Date.Collected, "%m/%d/%Y") == eventDate
                & Transect == transect) %>%
          arrange(Sample.Number) %>%
          pull(Sample.Number)

        site_count <- length(sites)
        for (k in 1:site_count) {
          site <- sites[k]
          eventID <- paste0(transect, "_", site)

          site_letter <- substr(site, nchar(site), nchar(site))

          species_count <- length(species)
          for (m in 1:species_count) {
            specie <- species[m]

            print(paste0("Transect(", i, "/", transect_count,
                        "); Date.Collected(", j, "/", collection_date_count,
                        "); Site(", k, "/", site_count,
                        "); Species(", m, "/", species_count, ")"))

            organism <- mycopins.organisms %>%
              filter(organism == specie) %>%
              slice(1)

            if (nrow(organism) == 1) { # exists
              occurrenceID <- paste0(eventID, ":", organism$gbif.usage_key)

              occurrence_match <- mycopins %>%
                filter(Transect == transect
                      & as.Date(Date.Collected, "%m/%d/%Y") == eventDate
                      & substrRight(Sample.Number, 1) == site_letter) %>%
                slice(1)

              occurrence_record <- create_gbif_occurrence_record(eventID, occurrenceID,
                                                            basisOfRecord, organismQuantityType, sampleSizeUnit,
                                                            site_letter, specie,
                                                            organism, occurrence_match)
              gbif_occurrences <- rbind(gbif_occurrences,
                                        as.data.frame(occurrence_record, stringsAsFactors = FALSE))
            } else {
              # Do nothing.
            }
          }
        }
      }
    }
  })

  print(execution_time)

  return(gbif_occurrences)
}

# > print(execution_time)
# user  system elapsed
# 300.69    5.06  841.31