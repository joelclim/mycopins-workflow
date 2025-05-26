if (!require("dplyr")) install.packages("dplyr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("readr")) install.packages("readr")

library(dplyr)
library(jsonlite)
library(readr)

#
# Requires get_transects(), get_habitat(transect), get_latitude(transect), get_longitude(transect)
#          get_wood_type(site_letter), get_wood_texture(site_letter), get_match_source(source)
#

create_gbif_occurrence_dynamic_properties <- function(site_letter, sample_cluster, organism) {
  dynamic_properties_df <- data.frame(woodType = get_wood_type(site_letter),
                                      woodTexture = get_wood_texture(site_letter))

  if (nrow(sample_cluster) > 0) {
    match_df <- sample_cluster[, c("source", "percent_identity",
                                  "evalue", "bit_score",
                                  "accession", "taxid",
                                  "Count")]

    dynamic_properties_df$match <- list(match_df)
  }

  if (nrow(organism) > 0 &
    'fungal_traits.primary_lifestyle'%in% colnames(organism)) {
    fungal_traits_df <- data.frame(
      primaryLifestyle = organism$fungal_traits.primary_lifestyle,
      secondaryLifestyle = organism$fungal_traits.secondary_lifestyle
    )
    dynamic_properties_df$fungalTraits <- fungal_traits_df
  }

  return(dynamic_properties_df)
}

create_gbif_occurrence_record <- function(eventID, occurrenceID,
                                    basisOfRecord, organismQuantityType,
                                    site_letter, specie, sample_cluster, organism, occurrence_match) {
  get_data_generalizations <- function(organism) {
    if (organism$source == "BLAST") {
      return(paste0("Identified via BLASTn circa 2024-05-01 with ",
                  organism$bit_score, " bit score, ",
                  organism$evalue, " evalue, and ",
                  format(round(organism$percent_identity, 2), nsmall = 2), "% identity."))
    }
    return("Identified via UNITE Fungi 9.0 (2023-07-18) (https://doi.org/10.15156/BIO/2938068)")
  }

  get_taxon_id <- function(organism) {
    return(paste0("https://www.gbif.org/species/", organism$gbif.usage_key))
  }

  get_taxon_concept_id <- function(organism) {
    if (organism$source == "BLAST") {
      return(paste0("NCBI:tx", organism$taxid))
    }

    # UNITE
    tax_id <- strsplit(organism$taxid, split="\\|")
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

  dynamic_properties_df <- create_gbif_occurrence_dynamic_properties(site_letter, sample_cluster, organism)

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
    dynamicProperties = toJSON(dynamic_properties_df)
  )

  return(occurrence_record)
}

#
# Main function: generate_gbif_occurrence
#
generate_gbif_event_occurrence <- function(configuration) {
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }

  location <- configuration["location"]
  basisOfRecord <- configuration["basisOfRecord"]
  organismQuantityType <- configuration["organismQuantityType"]

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

  sample_cluster_file <- configuration["mycopins_sample_cluster_file"]
  mycopins.sample_cluster <- read_csv(sample_cluster_file, show_col_types = FALSE,
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

        samples <- mycopins.environment %>%
          filter(as.Date(Date.Collected, "%m/%d/%Y") == eventDate
                & Transect == transect) %>%
          arrange(Sample.Number) %>%
          pull(Sample.Number)

        # site = sample
        sample_count <- length(samples)
        for (k in 1:sample_count) {
          sample <- samples[k]
          eventID <- paste0(transect, "_", sample)

          # For wood type and wood texture classification
          site_letter <- substr(sample, nchar(sample), nchar(sample))

          species_count <- length(species)
          for (m in 1:species_count) {
            specie <- species[m]

            if (specie %in% c('mock')) {
              next
            }
            # if (specie %in% c('Fungi', 'Fungus', 'mock')) {
            #   next
            # }

            print(paste0("Transect(", i, "/", transect_count,
                        "); Date.Collected(", j, "/", collection_date_count,
                        "); Sample(", k, "/", sample_count,
                        "); Species(", m, "/", species_count, ")"))


            sample_cluster <- mycopins.sample_cluster %>%
              filter(Transect == transect &
                     Sample.Number == sample &
                     gbif_accepted_name == specie & gbif.kingdom == 'Fungi') %>%
              arrange(-Count)

            organism <- NULL
            if (nrow(sample_cluster) == 0) {
              # absent: first instance of the specie regardless of the transect
              organism <- mycopins.sample_cluster %>%
                            filter(gbif_accepted_name == specie & gbif.kingdom == 'Fungi') %>%
                            dplyr::slice(1)
            } else {
              # present: The most count of the selected sample_cluster records
              organism <- sample_cluster[1, ]
            }

            # absent, gbif.kingdom != 'Fungi'
            if (nrow(organism) != 1) {
              next
            }

            occurrenceID <- paste0(eventID, ":", organism$gbif.usage_key)
            occurrence_match <- mycopins %>%
              filter(Transect == transect
                    & Sample.Number == sample
                    & as.Date(Date.Collected, "%m/%d/%Y") == eventDate)

            occurrence_record <- create_gbif_occurrence_record(eventID, occurrenceID,
                                                          basisOfRecord, organismQuantityType,
                                                          site_letter, specie, sample_cluster,
                                                          organism, occurrence_match)
            gbif_occurrences <- rbind(gbif_occurrences,
                                      as.data.frame(occurrence_record, stringsAsFactors = FALSE))
          }
        }
      }
    }
  })

  print(execution_time)

  return(gbif_occurrences)
}


generate_gbif_event_occurrence_no_kingdom <- function(configuration) {
  ipt_gbif_event_occurrence_all_fungi_file <- configuration["ipt_gbif_event_occurrence_all_fungi_file"]
  gbif.event_occurrence <- read_csv(ipt_gbif_event_occurrence_all_fungi_file, show_col_types = FALSE,
                                  locale = locale(encoding = "UTF-8"))



}

# > print(execution_time)
# user  system elapsed
# 300.69    5.06  841.31

#    user  system elapsed
# 5243.58   66.36 5628.10