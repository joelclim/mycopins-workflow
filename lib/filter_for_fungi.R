filter_for_fungi <- function(counts_df,
                             organism_df,
                             sep_col_name = "_Splitter_") {
  fungi_taxon_df <- organism_df[nchar(organism_df$gbif.genus) > 0 &
                                organism_df$gbif.kingdom == "Fungi",]

  start_col <- which(names(counts_df) == sep_col_name)

  wood.env <- names(counts_df[, 1:start_col])

  # For each record, get the value in acceptedNameUsage if exists.
  # Otherwise, use the originalNameUsage.
  fungi_names <- unique(sort(ifelse(!is.na(fungi_taxon_df$gbif.accepted_name_usage),
                                    fungi_taxon_df$gbif.accepted_name_usage,
                                    fungi_taxon_df$gbif.original_name_usage)))

  # TODO: QC
  # Exclude: fungi with names like "Fungus", "Fungi", or mock
  fungi_names <- setdiff(fungi_names, c("alga",
                                        "fungi", "fungus",
                                        "Fungi", "Fungus",
                                        "mock", "uncultured organism"))

  # TODO: QC
  # Exclude: fungi with bit_score less than 200
  less_than_200_df <- subset(organism_df, match.bit_score < 200)
  less_than_200_fungi_names <- unique(sort(ifelse(!is.na(less_than_200_df$gbif.accepted_name_usage),
                                    less_than_200_df$gbif.accepted_name_usage,
                                    less_than_200_df$gbif.original_name_usage)))

  fungi_names <- setdiff(fungi_names, less_than_200_fungi_names)

  fungi_only_df <- counts_df[names(counts_df) %in% c(wood.env, fungi_names)]

  # Remove rows whose sum of counts is zero
  num_counts_df <- fungi_only_df[, (start_col+1):ncol(fungi_only_df)]
  fungi_only_df <- fungi_only_df[rowSums(num_counts_df) != 0, ]

  return(fungi_only_df)
}
