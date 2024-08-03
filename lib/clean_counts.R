if (!require("readr")) install.packages("readr")

library(readr)

clean_counts <- function(scata_counts_file,
                         env_file,
                         dataset_name,
                         cleaned_counts_qc_file,
                         minimum_count = 3,
                         sep_col_name = "_Splitter_") {

  # Keep tags where the forward primer matches with the reverse primer.
  df <- clean_counts.valid_tags(scata_counts_file, dataset_name)
  counts <- df
  counts.tags <- counts$Tag

  # Load the environment file
  env_df <- read.csv(env_file)
  env.tags <- env_df$Tag

  # Tags in env but not in counts
  env.tags_not_in_counts <- setdiff(env.tags, counts.tags)

  # Tags in counts but not in env
  counts.tags_not_in_env <- setdiff(counts.tags, env.tags)
  counts_not_in_env_df <- counts[counts$Tag %in% counts.tags_not_in_env, ]

  # Keep tags that exists in both env and counts
  shared_tags <- intersect(env.tags, counts.tags)
  env_df <- env_df[env_df$Tag %in% shared_tags, ]
  counts <- counts[counts$Tag %in% shared_tags, ]
  counts$Sample.Number <- env_df$Sample.Number

  cluster_counts <- ncol(counts) - 2 # Exclude Tag and Sample.Number

  counts.cols_with_pos <- character()
  if (any(counts$Sample.Number == "POS")) {
  # Remove clusters with positive values greater than minimum count
    counts.cols_with_pos <- clean_counts.cols_with_control(counts, "POS", minimum_count)
    counts <- counts[, !(names(counts) %in% counts.cols_with_pos)]
    counts <- counts[counts$Sample.Number != "POS", ]
    env_df <- env_df[env_df$Sample.Number != "POS", ]
  }

  counts.cols_with_neg <- character()
  if (any(counts$Sample.Number == "NEG")) {
    # Remove clusters with negative values greater than minimum count
    counts.cols_with_neg <- clean_counts.cols_with_control(counts, "NEG", minimum_count)
    counts <- counts[, !(names(counts) %in% counts.cols_with_neg)]
    counts <- counts[counts$Sample.Number != "NEG", ]
    env_df <- env_df[env_df$Sample.Number != "NEG", ]
  }

  # Take the community data from counts
  community_df <- counts[, !(names(counts) %in% c("Tag", "Sample.Number")) ]

  # Exclude small clusters
  small_clusters <- clean_counts.get_small_clusters(community_df, minimum_count)
  community_df <- community_df[, !(names(community_df) %in% small_clusters)]

  qc_lines <- c(
    paste0("Final number of env tags (A - C - {POS, NEG}) = ", nrow(env_df)),
    paste0("Final number of clusters (E - F - G - H) = ", ncol(community_df)),
    paste(rep("=", 80), collapse = ""),
    paste0("[A] Number of env tags: ", length(env.tags)),
    paste0("[B] Number of counts tags: ", length(counts.tags)),
    paste0("[C] Number of env tags not in counts: ", length(env.tags_not_in_counts)),
    paste0("[D] Number of counts tags not in env: ", length(counts.tags_not_in_env)),
    paste0("[E] Number of clusters: ", cluster_counts),
    paste0("[F] Number of clusters with pos control: ", length(counts.cols_with_pos)),
    paste0("[G] Number of clusters with neg control: ", length(counts.cols_with_neg)),
    paste0("[H] Number of clusters with total counts less than ",
           minimum_count, ": ", length(small_clusters)),
    paste(rep("=", 80), collapse = ""),
    "Tags in env not in counts. Removed from final env tags.",
    ifelse(0 == length(env.tags_not_in_counts), "None",
                                                paste(env.tags_not_in_counts, collapse = ", ")),
    paste(rep("=", 80), collapse = ""),
    "Tags in counts not in env. Excluded in final env tags.",
    ifelse(0 == length(counts.tags_not_in_env), "None",
                                                paste(counts.tags_not_in_env, collapse = ", ")),
    paste(rep("=", 80), collapse = ""),
    "Clusters with POS control",
    ifelse(0 == length(counts.cols_with_pos), "None",
                                              paste(counts.cols_with_pos, collapse = ", ")),
    paste(rep("=", 80), collapse = ""),
    "Clusters with NEG control",
    ifelse(0 == length(counts.cols_with_neg), "None",
                                              paste(counts.cols_with_neg, collapse = ", ")),
    paste(rep("=", 80), collapse = ""),
    paste("Clusters with total counts less than", minimum_count),
    ifelse(0 == length(small_clusters), "None",
                                         paste(small_clusters, collapse = ", "))
  )
  writeLines(qc_lines, cleaned_counts_qc_file)

  sep_df <- data.frame(matrix(ncol = 1, nrow = 1))
  colnames(sep_df) <- sep_col_name
  sep_df[, 1] <- "X"

  return(cbind(env_df, sep_df, community_df))
}


clean_counts.valid_tags <- function(in_counts_file, dataset_name) {
  # Load the counts file
  df <- read.csv(in_counts_file, sep = ";")

  # Filter the dataset based on the pattern
  pattern <- paste(dataset_name, "tag([^\\s]+)([FR])_tag([^\\s]+)([FR])")
  df_tags <- df[grepl(pattern, df$Tag), ]

  # Extract group 2 and group 3 from the matched rows
  # ([^\\s]+)_([^\\s]+) groups 2 and 3
  df_tags$TagNum1 <- as.numeric(sub(pattern, "\\1", df_tags$Tag))
  df_tags$TagNum1FR <- sub(pattern, "\\2", df_tags$Tag)
  df_tags$TagNum2 <- as.numeric(sub(pattern, "\\3", df_tags$Tag))
  df_tags$TagNum2FR <- sub(pattern, "\\4", df_tags$Tag)

  # Sort by tag numbers
  sorted_df_tags <- df_tags[order(df_tags$TagNum1,
                                  df_tags$TagNum1FR,
                                  df_tags$TagNum2,
                                  df_tags$TagNum2FR), ]

  # Remove entries where the forward tag num does not match
  # with the reverse tag num.
  df_valid_tags <- sorted_df_tags[
    sorted_df_tags$TagNum1 == sorted_df_tags$TagNum2, ]

  # Remove columns we no longer need
  df_valid_tags <- df_valid_tags[, !(names(df_valid_tags) %in%
                                       c("Tag", "TagNum2",
                                         "TagNum1FR", "TagNum2FR"))]

  # Move TagNum1 column as the first column.
  df_valid_tags <- df_valid_tags[, c("TagNum1",
                                     setdiff(names(df_valid_tags), "TagNum1"))]

  # Rename TagNum1 as Tag
  names(df_valid_tags)[names(df_valid_tags) == "TagNum1"] <- "Tag"

  return(df_valid_tags)
}

clean_counts.get_small_clusters <- function(community_df, minimum_count = 3) {
  # Calculate the total value of each column
  col_sums <- colSums(community_df)

  # Find the names of columns where the total value is less than the minimum count
  small_clusters <- names(col_sums[col_sums < minimum_count])

  return(small_clusters)
}

clean_counts.cols_with_control <- function(counts, control, minimum_count = 3) {
  control_rows <- counts[counts$Sample.Number == control,
                            !(colnames(counts) %in% c("Tag", "Sample.Number"))]
  col_names <- colnames(control_rows)
  col_names_with_control <- col_names[apply(control_rows, 2,
                                            function(x) max(x) >= minimum_count)]

  return(col_names_with_control)
}
