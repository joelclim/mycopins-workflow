filter_for_fungi <- function(counts_df,
                             gbif_taxon_df,
                             sep_col_name = "_Splitter_") {
  fungi_taxon_df <- gbif_taxon_df[nchar(gbif_taxon_df$genus) > 0,]

  start_col <- which(names(counts_df) == sep_col_name)

  wood.env <- names(counts_df[, 1:start_col])

  # For each record, get the value in acceptedNameUsage if exists.
  # Otherwise, use the originalNameUsage.
  fungi_names <- unique(sort(ifelse(!is.na(fungi_taxon_df$acceptedNameUsage),
                                    fungi_taxon_df$acceptedNameUsage,
                                    fungi_taxon_df$originalNameUsage)))

  # Exclude: "fungus" in ('uncultured fungi', 'fungal sp.')
  fungi_only_df <- counts_df[names(counts_df) %in% c(wood.env, fungi_names)]
  fungi_only_df <- fungi_only_df[names(fungi_only_df) != "Fungus"]

  # Remove rows whose sum of counts is zero
  num_counts_df <- fungi_only_df[, (start_col+1):ncol(fungi_only_df)]
  fungi_only_df <- fungi_only_df[rowSums(num_counts_df) != 0, ]

  return(fungi_only_df)
}