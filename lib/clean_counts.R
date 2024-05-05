#
# A script that creates a cleaned version of the SCATA all_tag_by_cluster_counts counts file.
# It also creates a table of cluster ids with a pipe-separated list of tags present in each cluster.
#


################################################################################
## Create a CSV file a cleaned version of the SCATA all_tag_by_cluster_counts.csv
## result file.Requires the count file and the label file.
##
## The cleaning steps includes:
## 1. Remove tag records that do not have matching forward and reverse tags
## 2. Keep clusters that contain more than 2 total values.
## 3. Remove clusters that have more than 2 positive control values.
## 4. Remove clusters that have more than 2 negative control values.
## 5. Combine labels and cluster data as a cleaned data frame.
## 6. Save the cleaned data frame as a CSV file.
## NOTE: Please provide the in_features_file name and in_counts_file name
################################################################################
clean_counts <- function(in_counts_file,
                         in_features_file,
                         dataset_name,
                         minimum_count=3,
                         sep_col_name = "_Splitter_") {
  # Load the counts file
  df <- read.csv(in_counts_file, sep = ";")

  # Load the labels file
  features_df <- read.csv(in_features_file)

  # Filter the dataset based on the pattern
  pattern <- paste(dataset_name, "tag([^\\s]+)([FR])_tag([^\\s]+)([FR])")
  df_tags <- df[grepl(pattern, df$Tag), ]

  # Extract group 2 and group 3 from the matched rows
  # ([^\\s]+)_([^\\s]+) groups 2 and 3
  df_tags$TagNum1 <- sub(pattern, "\\1", df_tags$Tag)
  df_tags$TagNum1FR <- sub(pattern, "\\2", df_tags$Tag)
  df_tags$TagNum2 <- sub(pattern, "\\3", df_tags$Tag)
  df_tags$TagNum2FR <- sub(pattern, "\\4", df_tags$Tag)

  sorted_df_tags <- df_tags[order(df_tags$TagNum1,
                                  df_tags$TagNum1FR,
                                  df_tags$TagNum2,
                                  df_tags$TagNum2FR), ]

  # Remove invalid tags
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

  # Split Tag column and cluster columns
  tags_df <- df_valid_tags[, c("Tag") ]
  clusters_df <- df_valid_tags[, !(names(df_valid_tags) %in% c("Tag")) ]

  # Calculate the total value of each cluster column
  column_totals <- colSums(clusters_df)
  # Find the names of columns where the total value is greater than 2
  selected_columns <- names(column_totals[column_totals >= minimum_count])
  # Subset the clusters dataframe to include only the selected columns
  clusters_df <- df_valid_tags[, selected_columns]
  # Name each row by the tag names
  rownames(clusters_df) <- tags_df

  discriminator_df <- data.frame(matrix(ncol = 1, nrow = 1))
  colnames(discriminator_df) <- sep_col_name
  discriminator_df[, 1] <- "X"

  # Combine tags, labels, clusters
  combined_df <- cbind(tags_df, features_df,
                       discriminator_df, clusters_df)

  # Determine if there are row tags that do not match with the label tags
  # NOTE: qc_match_df should have 0 observations. Otherwise, DO NOT CONTINUE.
  # TODO: Use qc_match_df for validation
  qc_match_df <- combined_df[combined_df$tags_df != combined_df$Tag, ]

  clean_df <- combined_df[, !(names(combined_df) %in% c("tags_df"))]

  # Remove columns with positive and negative control values greater than 2
  pos_rows <- clean_df[clean_df$Sample.Number == "POS",
                       !(names(clean_df) %in%
                           c(names(features_df), sep_col_name))]
  columns_with_pos <- colnames(pos_rows)[apply(pos_rows, 2, function(x) max(x) >= minimum_count)]
  clean_df <- clean_df[, !(names(clean_df) %in% columns_with_pos)]

  neg_rows <- clean_df[clean_df$Sample.Number == "NEG",
                       !(names(clean_df) %in%
                           c(names(features_df), sep_col_name))]
  columns_with_neg <- colnames(neg_rows)[apply(neg_rows, 2, function(x) max(x) >= minimum_count)]
  clean_df <- clean_df[, !(names(clean_df) %in% columns_with_neg)]

  # Remove the positive and negative control value rows
  clean_df <- clean_df[clean_df$Sample.Number != "POS", ]
  clean_df <- clean_df[clean_df$Sample.Number != "NEG", ]

  return(clean_df)
}
