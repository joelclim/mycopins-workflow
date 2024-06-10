if (!require("dplyr")) install.packages("dplyr")
if (!require("rgbif")) install.packages("rgbif")

library("dplyr")
library("rgbif")

get_gbif_taxon <- function(organisms, complete_clusters_df, fungal_traits_df) {
  get_gbif_taxon_by_organism <- function(organism, complete_clusters_df, fungal_traits_df) {
    response <- name_backbone_verbose(name=organism, kingdom="fungi")
    taxon <- response$data
    if (taxon$matchType == "NONE" | organism == "mock") {
      return(c(
        searchSource = "",
        searchPercentIdentity = "",
        searchEValue = "",
        searchBitScore = "",
        searchAccession = "",
        searchTaxId = "",
        matchType = taxon$matchType,
        scientificName = "",
        acceptedNameUsage = "",
        originalNameUsage = "",
        kingdom = "",
        phylum = "",
        class = "",
        order = "",
        family = "",
        genus = "",
        taxonKey = "",
        status = "",
        scientificNameAuthorship = "",
        fungalTraitsPrimaryLifestyle = "",
        fungalTraitsSecondaryLifestyle = ""
      ))
    }

    # Match search results: UNITE or BLAST
    search_matches <- complete_clusters_df %>%
      filter(normalized_name == organism) %>%
      arrange(desc(percent_identity), evalue)

    best_search_match <- NULL
    if (nrow(search_matches) > 0) {
      best_search_match <- search_matches[1,]
    }

    # Fungal traits of the genus
    genus_traits <- NULL
    if ("genus" %in% names(taxon)) {
      genus_traits <- fungal_traits_df[fungal_traits_df$GENUS == taxon$genus,]
    }

    gbif_taxon <- c(
      searchSource = ifelse(nrow(best_search_match) == 1, best_search_match$source, ""),
      searchPercentIdentity = ifelse(nrow(best_search_match) == 1, best_search_match$percent_identity, ""),
      searchEValue = ifelse(nrow(best_search_match) == 1,  best_search_match$evalue, ""),
      searchBitScore = ifelse(nrow(best_search_match) == 1, best_search_match$bit_score, ""),
      searchAccession = ifelse(nrow(best_search_match) == 1, best_search_match$accession, ""),
      searchTaxId = ifelse(nrow(best_search_match) == 1,  best_search_match$taxid, ""),
      matchType = ifelse("matchType" %in% names(taxon),
                         taxon$matchType, ""),
      scientificName = ifelse("scientificName" %in% names(taxon),
                              taxon$scientificName, ""),
      acceptedNameUsage = ifelse("species" %in% names(taxon),
                                 taxon$species, ""),
      originalNameUsage = ifelse("canonicalName" %in% names(taxon),
                                 taxon$canonicalName, ""),
      kingdom = ifelse("kingdom" %in% names(taxon), taxon$kingdom, ""),
      phylum = ifelse("phylum" %in% names(taxon), taxon$phylum, ""),
      class = ifelse("class" %in% names(taxon), taxon$class, ""),
      order = ifelse("order" %in% names(taxon), taxon$order, ""),
      family = ifelse("family" %in% names(taxon), taxon$family, ""),
      genus = ifelse("genus" %in% names(taxon), taxon$genus, ""),
      taxonKey = ifelse("rank" %in% names(taxon), taxon$rank, ""),
      status = ifelse("status" %in% names(taxon), taxon$status, ""),
      scientificNameAuthorship = ifelse(
        "scientificName" %in% names(taxon) &
        "canonicalName" %in% names(taxon),
        trimws(substring(
          taxon$scientificName,
          nchar(taxon$canonicalName)+1,
          nchar(taxon$scientificName))), ""),
      fungalTraitsPrimaryLifestyle = ifelse(nrow(genus_traits) == 1, genus_traits$primary_lifestyle, ""),
      fungalTraitsSecondaryLifestyle = ifelse(nrow(genus_traits) == 1, genus_traits$Secondary_lifestyle, "")
    )

    return(gbif_taxon)
  }

  organism_gbif_df <- data.frame(organism = organisms)

  searchSource <- character(nrow(organism_gbif_df))
  searchPercentIdentity <- numeric(nrow(organism_gbif_df))
  searchEValue <- numeric(nrow(organism_gbif_df))
  searchBitScore <- numeric(nrow(organism_gbif_df))
  searchAccession <- character(nrow(organism_gbif_df))
  searchTaxId <- character(nrow(organism_gbif_df))
  matchType <- character(nrow(organism_gbif_df))
  scientificName <- character(nrow(organism_gbif_df))
  acceptedNameUsage <- character(nrow(organism_gbif_df))
  originalNameUsage <- character(nrow(organism_gbif_df))
  kingdom <- character(nrow(organism_gbif_df))
  phylum <- character(nrow(organism_gbif_df))
  class <- character(nrow(organism_gbif_df))
  order <- character(nrow(organism_gbif_df))
  family <- character(nrow(organism_gbif_df))
  genus <- character(nrow(organism_gbif_df))
  taxonKey <- character(nrow(organism_gbif_df))
  status <- character(nrow(organism_gbif_df))
  scientificNameAuthorship <- character(nrow(organism_gbif_df))
  fungalTraitsPrimaryLifestyle <- character(nrow(organism_gbif_df))
  fungalTraitsSecondaryLifestyle <- character(nrow(organism_gbif_df))

  for (i in 1:nrow(organism_gbif_df)) {
    taxon <- get_gbif_taxon_by_organism(organism_gbif_df[i, ], complete_clusters_df, fungal_traits_df)

    searchSource[i] <- taxon['searchSource']
    searchPercentIdentity[i] <- taxon['searchPercentIdentity']
    searchEValue[i] <- taxon['searchEValue']
    searchBitScore[i] <- taxon['searchBitScore']
    searchAccession[i] <- taxon['searchAccession']
    searchTaxId[i] <- taxon['searchTaxId']
    matchType[i] <- taxon['matchType']
    scientificName[i] <- taxon['scientificName']
    acceptedNameUsage[i] <- taxon['acceptedNameUsage']
    originalNameUsage[i] <- taxon['originalNameUsage']
    kingdom[i] <- taxon['kingdom']
    phylum[i] <- taxon['phylum']
    class[i] <- taxon['class']
    order[i] <- taxon['order']
    family[i] <- taxon['family']
    genus[i] <- taxon['genus']
    taxonKey[i] <- taxon['taxonKey']
    status[i] <- taxon['status']
    scientificNameAuthorship[i] <- taxon['scientificNameAuthorship']
    fungalTraitsPrimaryLifestyle[i] <- taxon['fungalTraitsPrimaryLifestyle']
    fungalTraitsSecondaryLifestyle[i] <- taxon['fungalTraitsSecondaryLifestyle']
  }

  organism_gbif_df$searchSource <- searchSource
  organism_gbif_df$searchPercentIdentity <- searchPercentIdentity
  organism_gbif_df$searchEValue <- searchEValue
  organism_gbif_df$searchBitScore <- searchBitScore
  organism_gbif_df$searchAccession <- searchAccession
  organism_gbif_df$searchTaxId <- searchTaxId
  organism_gbif_df$matchType <- matchType
  organism_gbif_df$scientificName <- scientificName
  organism_gbif_df$acceptedNameUsage <- acceptedNameUsage
  organism_gbif_df$originalNameUsage <- originalNameUsage
  organism_gbif_df$kingdom <- kingdom
  organism_gbif_df$phylum <- phylum
  organism_gbif_df$class <- class
  organism_gbif_df$order <- order
  organism_gbif_df$family <- family
  organism_gbif_df$genus <- genus
  organism_gbif_df$taxonKey <- taxonKey
  organism_gbif_df$status <- status
  organism_gbif_df$scientificNameAuthorship <- scientificNameAuthorship
  organism_gbif_df$fungalTraitsPrimaryLifestyle <- fungalTraitsPrimaryLifestyle
  organism_gbif_df$fungalTraitsSecondaryLifestyle <- fungalTraitsSecondaryLifestyle

  return(organism_gbif_df)
}
