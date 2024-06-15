if (!require("dplyr")) install.packages("dplyr")
if (!require("rgbif")) install.packages("rgbif")

library("dplyr")
library("rgbif")

# organisms is a list of gbif_accepted_names
create_organism_dataset <- function(organisms, complete_clusters_df, fungal_traits_df) {
  create_organism_record <- function(organism, complete_clusters_df, fungal_traits_df) {
    response <- name_backbone_verbose(name=organism, kingdom="fungi")
    taxon <- response$data

    if (taxon$matchType == "NONE" | organism == "mock") {
      return(c(
        match.source = "",
        match.percent_identity = "",
        match.evalue = "",
        match.bit_score = "",
        match.accession = "",
        match.tax_id = "",
        gbif.match_type = taxon$matchType,
        gbif.usage_key = "",
        gbif.scientific_name = "",
        gbif.accepted_name_usage = "",
        gbif.original_name_usage = "",
        gbif.kingdom = "",
        gbif.phylum = "",
        gbif.class = "",
        gbif.order = "",
        gbif.family = "",
        gbif.genus = "",
        gbif.taxon_rank = "",
        gbif.status = "",
        gbif.scientific_name_authorship = "",
        fungal_traits.primary_lifestyle = "",
        fungal_traits.secondary_lifestyle = ""
      ))
    }

    # Match search results: UNITE or BLAST
    search_matches <- complete_clusters_df %>%
      filter(gbif_accepted_name == organism) %>%
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

    organism_record <- c(
      match.source = ifelse(nrow(best_search_match) == 1,
                             best_search_match$source, ""),
      match.percent_identity = ifelse(nrow(best_search_match) == 1,
                                      best_search_match$percent_identity, ""),
      match.evalue = ifelse(nrow(best_search_match) == 1,
                             best_search_match$evalue, ""),
      match.bit_score = ifelse(nrow(best_search_match) == 1,
                               best_search_match$bit_score, ""),
      match.accession = ifelse(nrow(best_search_match) == 1,
                                best_search_match$accession, ""),
      match.tax_id = ifelse(nrow(best_search_match) == 1,
                             best_search_match$taxid, ""),
      gbif.match_type = ifelse("matchType" %in% names(taxon),
                         taxon$matchType, ""),
      gbif.usage_key = ifelse("usageKey" %in% names(taxon),
                         taxon$usageKey, ""),
      gbif.scientific_name = ifelse("scientificName" %in% names(taxon),
                              taxon$scientificName, ""),
      gbif.accepted_name_usage = ifelse("species" %in% names(taxon),
                                 taxon$species, ""),
      gbif.original_name_usage = ifelse("canonicalName" %in% names(taxon),
                                 taxon$canonicalName, ""),
      gbif.kingdom = ifelse("kingdom" %in% names(taxon), taxon$kingdom, ""),
      gbif.phylum = ifelse("phylum" %in% names(taxon), taxon$phylum, ""),
      gbif.class = ifelse("class" %in% names(taxon), taxon$class, ""),
      gbif.order = ifelse("order" %in% names(taxon), taxon$order, ""),
      gbif.family = ifelse("family" %in% names(taxon), taxon$family, ""),
      gbif.genus = ifelse("genus" %in% names(taxon), taxon$genus, ""),
      gbif.taxon_rank = ifelse("rank" %in% names(taxon), taxon$rank, ""),
      gbif.status = ifelse("status" %in% names(taxon), taxon$status, ""),
      gbif.scientific_name_authorship = ifelse(
        "scientificName" %in% names(taxon) &
        "canonicalName" %in% names(taxon),
        trimws(substring(
          taxon$scientificName,
          nchar(taxon$canonicalName)+1,
          nchar(taxon$scientificName))), ""),
      fungal_traits.primary_lifestyle = ifelse(nrow(genus_traits) == 1, genus_traits$primary_lifestyle, ""),
      fungal_traits.secondary_lifestyle = ifelse(nrow(genus_traits) == 1, genus_traits$Secondary_lifestyle, "")
    )

    return(organism_record)
  }

  organism_df <- data.frame(organism = organisms)

  match.source <- character(nrow(organism_df))
  match.percent_identity <- numeric(nrow(organism_df))
  match.evalue <- numeric(nrow(organism_df))
  match.bit_score <- numeric(nrow(organism_df))
  match.accession <- character(nrow(organism_df))
  match.tax_id <- character(nrow(organism_df))
  gbif.match_type <- character(nrow(organism_df))
  gbif.usage_key <- character(nrow(organism_df))
  gbif.scientific_name <- character(nrow(organism_df))
  gbif.accepted_name_usage <- character(nrow(organism_df))
  gbif.original_name_usage <- character(nrow(organism_df))
  gbif.kingdom <- character(nrow(organism_df))
  gbif.phylum <- character(nrow(organism_df))
  gbif.class <- character(nrow(organism_df))
  gbif.order <- character(nrow(organism_df))
  gbif.family <- character(nrow(organism_df))
  gbif.genus <- character(nrow(organism_df))
  gbif.taxon_rank <- character(nrow(organism_df))
  gbif.status <- character(nrow(organism_df))
  gbif.scientific_name_authorship <- character(nrow(organism_df))
  fungal_traits.primary_lifestyle <- character(nrow(organism_df))
  fungal_traits.secondary_lifestyle <- character(nrow(organism_df))

  for (i in 1:nrow(organism_df)) {
    taxon <- create_organism_record(organism_df[i, ], complete_clusters_df, fungal_traits_df)

    match.source[i] <- taxon['match.source']
    match.percent_identity[i] <- taxon['match.percent_identity']
    match.evalue[i] <- taxon['match.evalue']
    match.bit_score[i] <- taxon['match.bit_score']
    match.accession[i] <- taxon['match.accession']
    match.tax_id[i] <- taxon['match.tax_id']
    gbif.match_type[i] <- taxon['gbif.match_type']
    gbif.usage_key[i] <- taxon['gbif.usage_key']
    gbif.scientific_name[i] <- taxon['gbif.scientific_name']
    gbif.accepted_name_usage[i] <- taxon['gbif.accepted_name_usage']
    gbif.original_name_usage[i] <- taxon['gbif.original_name_usage']
    gbif.kingdom[i] <- taxon['gbif.kingdom']
    gbif.phylum[i] <- taxon['gbif.phylum']
    gbif.class[i] <- taxon['gbif.class']
    gbif.order[i] <- taxon['gbif.order']
    gbif.family[i] <- taxon['gbif.family']
    gbif.genus[i] <- taxon['gbif.genus']
    gbif.taxon_rank[i] <- taxon['gbif.taxon_rank']
    gbif.status[i] <- taxon['gbif.status']
    gbif.scientific_name_authorship[i] <- taxon['gbif.scientific_name_authorship']
    fungal_traits.primary_lifestyle[i] <- taxon['fungal_traits.primary_lifestyle']
    fungal_traits.secondary_lifestyle[i] <- taxon['fungal_traits.secondary_lifestyle']
  }

  organism_df$match.source <- match.source
  organism_df$match.percent_identity <- match.percent_identity
  organism_df$match.evalue <- match.evalue
  organism_df$match.bit_score <- match.bit_score
  organism_df$match.accession <- match.accession
  organism_df$match.tax_id <- match.tax_id
  organism_df$gbif.match_type <- gbif.match_type
  organism_df$gbif.usage_key <- gbif.usage_key
  organism_df$gbif.scientific_name <- gbif.scientific_name
  organism_df$gbif.accepted_name_usage <- gbif.accepted_name_usage
  organism_df$gbif.original_name_usage <- gbif.original_name_usage
  organism_df$gbif.kingdom <- gbif.kingdom
  organism_df$gbif.phylum <- gbif.phylum
  organism_df$gbif.class <- gbif.class
  organism_df$gbif.order <- gbif.order
  organism_df$gbif.family <- gbif.family
  organism_df$gbif.genus <- gbif.genus
  organism_df$gbif.taxon_rank <- gbif.taxon_rank
  organism_df$gbif.status <- gbif.status
  organism_df$gbif.scientific_name_authorship <- gbif.scientific_name_authorship
  organism_df$fungal_traits.primary_lifestyle <- fungal_traits.primary_lifestyle
  organism_df$fungal_traits.secondary_lifestyle <- fungal_traits.secondary_lifestyle

  return(organism_df)
}
