filter_for_wood_saprotroph <- function(counts_df,
                                       gbif_taxon_df,
                                       wood_saprotroph_traits_df,
                                       sep_col_name = "_Splitter_") {
  fungi_taxon_df <- gbif_taxon_df[nchar(gbif_taxon_df$phylum) > 0,]
  wood_saprotroph_df <- fungi_taxon_df[fungi_taxon_df$genus
                                       %in% wood_saprotroph_traits_df$GENUS, ]

  start_col <- which(names(counts_df) == sep_col_name)

  wood.env <- names(counts_df[, 1:start_col])
  wood_saprotroph_names <- unique(sort(wood_saprotroph_df$organism))

  wood_saprotroph_only_df <- counts_df[names(counts_df)
                                       %in% c(wood.env, wood_saprotroph_names)]

  # Remove rows whose sum of counts is zero
  num_counts_df <- wood_saprotroph_only_df[, (start_col+1):ncol(wood_saprotroph_only_df)]
  wood_saprotroph_only_df <- wood_saprotroph_only_df[rowSums(num_counts_df) != 0, ]

  return(wood_saprotroph_only_df)
}