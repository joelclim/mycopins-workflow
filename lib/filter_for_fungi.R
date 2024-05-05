filter_for_fungi <- function(counts_df, 
                             gbif_taxon_df, 
                             sep_col_name = "_Splitter_") {
  fungi_df <- gbif_taxon_df[nchar(gbif_taxon_df$phylum) > 0,]
  
  start_col <- which(names(counts_df) == sep_col_name)
  
  features <- names(counts_df[, 1:start_col])
  
  # For each record, get the value in acceptedNameUsage if exists. 
  # Otherwise, get the value in organism (normalize name).
  fungi <- unique(sort(ifelse(nchar(fungi_df$acceptedNameUsage) > 0, 
                              fungi_df$acceptedNameUsage, fungi_df$originalNameUsage)
                       ))

  # fungus in ('uncultured fungi', 'fungal sp.')
  filtered_counts_df <- counts_df[names(counts_df) %in% c(features, fungi, "fungus")]
  
  return(filtered_counts_df)
}