apply_top_match_to_clusters <- function(top_match_df, cleaned_clusters_df) {
  names(cleaned_clusters_df)[names(cleaned_clusters_df) == "Cluster ID"] <- "cluster_id"

  # Merging data frames
  # Left-outer Join of the cleaned_clusters_df and top_match data frames.
  merge_df <- merge(cleaned_clusters_df, top_match_df,
                    by = "cluster_id", all.x = TRUE)

  # Determine the organism representing the cluster
  reference_df <- data.frame(t(apply(merge_df, 1, function(row) {
    # Use the top_match's sciname if exists
    if (!is.na(row["sciname"])) {
      return(c(row["sciname"],
               row["accession"],
               row["taxid"]))
    }

    # If not, use the cleaned_clusters_df's Reference.
    tokens <- unlist(strsplit(row["Reference"], split = "\\|"))

    # Specie's name of the organism is the first token of the Reference.
    sciname <- tokens[1]
    # Replace underscore with a space character.
    sciname <- gsub("_", " ", sciname)

    accession <- tokens[2]
    taxId <- paste0(tokens[3], "|", tokens[4])

    # Return the processed values
    return(c(sciname, accession, taxId))
  })))
  colnames(reference_df) <- c("sciname", "accession", "taxid")
  merge_df$Reference <- reference_df$sciname
  merge_df$normalized_name <- sapply(merge_df$Reference, normalize_name)
  merge_df$gbif_accepted_name <- sapply(merge_df$normalized_name,
                                       get_gbif_accepted_name)
  merge_df$accession <- reference_df$accession
  merge_df$taxid <- reference_df$taxid

  # Determine the whether SCATA or BLAST classified the cluster
  sources <- apply(merge_df, 1, function(row) {
    if (!is.na(row["sciname"])) {
      return("BLAST")
    }
    return("UNITE")
  })
  merge_df$source <- sources

  # Determine the cluster number for each cluster
  # Assumes cluster id is <job_id>_<cluster_id>.
  # For example: scata6390_1
  merge_df$cluster_num <- apply(merge_df, 1, function(row) {
    tokens <- strsplit(row["cluster_id"], "_")[[1]]
    return(as.numeric(tokens[length(tokens)]))
  })

  # Sort by cluster number
  merge_df <- merge_df[order(merge_df$cluster_num), ]

  # Keep all of the columns from cleaned_clusters_df plus the Organism column.
  keep_column_names <- names(cleaned_clusters_df)
  keep_column_names <- append(keep_column_names, "taxid", after = 4)
  keep_column_names <- append(keep_column_names, "accession", after = 4)
  keep_column_names <- append(keep_column_names, "source", after = 4)
  keep_column_names <- append(keep_column_names, "percent_identity", after = 4)
  keep_column_names <- append(keep_column_names, "gbif_accepted_name", after = 4)
  keep_column_names <- append(keep_column_names, "normalized_name", after = 4)
  merge_df <- merge_df[keep_column_names]

  # Drop all columns from top_match
  drop_column_names <- c(names(top_match_df), "X")
  drop_column_names <- drop_column_names[!drop_column_names %in%
                                           c("cluster_id",
                                             "percent_identity",
                                             "taxid",
                                             "accession")]
  merge_df <- merge_df[, !names(merge_df) %in% drop_column_names]

  # Percent Identity of Reference identified by SCATA is 100%S
  merge_df[["percent_identity"]][is.na(merge_df[["percent_identity"]])] <- 100

  # Restore the "Cluster ID" column name
  names(merge_df)[names(merge_df) == "cluster_id"] <- "Cluster.ID"

  return(merge_df)
}
