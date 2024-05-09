if (!require("readr")) install.packages("readr")

library(readr)

apply_match_to_counts <- function(counts_df,
                                  clusters_df,
                                  sep_col_name = "_Splitter_") {
  start_col <- which(names(counts_df) == sep_col_name)
  new_counts_df <- counts_df[, 1:start_col]

  names <- unique(sort(clusters_df$gbif_accepted_name))
  for (name in names) {
    clusters <- clusters_df$Cluster.ID[clusters_df$gbif_accepted_name == name]
    row_sum <- rowSums(counts_df[names(counts_df) %in% clusters])
    if (sum(row_sum) > 0) {
      new_counts_df[[name]] <- row_sum
    }
  }

  return(new_counts_df)
}