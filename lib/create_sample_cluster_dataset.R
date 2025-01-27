if (!require("rgbif")) install.packages("rgbif")
if (!require("tidyr")) install.packages("tidyr")

library("rgbif")
library("tidyr")

get_gbif_taxonomy_by_accepted_name <- function(gbif_accepted_name) {
  response <- name_backbone_verbose(name=gbif_accepted_name, kingdom="fungi")
  taxon <- response$data

  gbif_taxon_df <- data.frame(
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
        nchar(taxon$scientificName))), "")
  )

  return(gbif_taxon_df)
}

get_fungal_traits_by_genus <- function(genus, fungal_traits_df) {
  genus_traits <- NULL
  if (!is.na(genus)) {
    genus_traits <- fungal_traits_df[fungal_traits_df$GENUS == genus,]

    if (nrow(genus_traits) == 1) {
      return (c(
        genus_traits$primary_lifestyle,
        genus_traits$Secondary_lifestyle
      ))
    }
  }
  return (c('', ''))
}

create_sample_cluster_dataset <- function(cleaned_counts_df, complete_clusters_df, fungal_traits_df) {
  pivot_cluster_id_df <- cleaned_counts_df %>%
    pivot_longer(
      cols = starts_with("scata"),  # Select columns to pivot
      names_to = "Cluster.ID",      # New column for column names
      values_to = "Count"           # New column for values
    ) %>%
    filter(Count > 0) %>%           # Retain only rows where count > 0
    select(Transect, Sample.Number, Forward.Primer, Reverse.Primer, Cluster.ID, Count, )

  sample_cluster_df <- pivot_cluster_id_df %>%
    inner_join(complete_clusters_df, by='Cluster.ID')

  #
  # GBIF Taxonomy
  #
  gbif_taxon_df <- t(apply(sample_cluster_df, 1,
    function(row) {
      gbif_taxon <- get_gbif_taxonomy_by_accepted_name(row['gbif_accepted_name'])
      return(c(
        gbif_taxon$gbif.match_type,
        gbif_taxon$gbif.usage_key,
        gbif_taxon$gbif.scientific_name,
        gbif_taxon$gbif.accepted_name_usage,
        gbif_taxon$gbif.original_name_usage,
        gbif_taxon$gbif.kingdom,
        gbif_taxon$gbif.phylum,
        gbif_taxon$gbif.class,
        gbif_taxon$gbif.order,
        gbif_taxon$gbif.family,
        gbif_taxon$gbif.genus,
        gbif_taxon$gbif.taxon_rank,
        gbif_taxon$gbif.status,
        gbif_taxon$gbif.scientific_name_authorship
      ))
    }))

  # Convert to data frame and give appropriate column names
  gbif_taxon_df <- data.frame(gbif_taxon_df)
  colnames(gbif_taxon_df) <- c(
    'gbif.match_type',
    'gbif.usage_key',
    'gbif.scientific_name',
    'gbif.accepted_name_usage',
    'gbif.original_name_usage',
    'gbif.kingdom',
    'gbif.phylum',
    'gbif.class',
    'gbif.order',
    'gbif.family',
    'gbif.genus',
    'gbif.taxon_rank',
    'gbif.status',
    'gbif.scientific_name_authorship'
  )

  # Bind the new columns to the original data frame
  sample_cluster_df <- cbind(sample_cluster_df, gbif_taxon_df)

  #
  # Fungal Traits
  #
  genus_traits_df <- t(apply(sample_cluster_df, 1,
    function(row) {
      return (get_fungal_traits_by_genus(row['gbif.genus'], fungal_traits_df))
    }
  ))

  # Convert to data frame and give appropriate column names
  genus_traits_df <- data.frame(genus_traits_df)
  colnames(genus_traits_df) <- c(
    'fungal_traits.primary_lifestyle',
    'fungal_traits.secondary_lifestyle'
  )

  # Bind the new columns to the original data frame
  sample_cluster_df <- cbind(sample_cluster_df, genus_traits_df)

  return(sample_cluster_df)
}
