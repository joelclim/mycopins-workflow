merge_add_counts <- function(top_df, bot_df,
                             sep_col_name = "_Splitter_") {
  start_col <- which(names(top_df) == sep_col_name)

  merged_df <- rbind(top_df[, 1:start_col], bot_df[,1:start_col])

  top_column_names <- unique(sort(c(names(
    top_df[, (start_col+1):ncol(top_df)]))))
  bot_column_names <- unique(sort(c(names(
    bot_df[, (start_col+1):ncol(bot_df)]))))
  col_names <- unique(sort(c(top_column_names, bot_column_names)))

  for (col_name in col_names) {
    top_values <- rep(0, nrow(top_df))
    bot_values <- rep(0, nrow(bot_df))

    if (col_name %in% top_column_names) {
      top_values <- top_df[, col_name]
    }
    if (col_name %in% bot_column_names) {
      bot_values <- bot_df[, col_name]
    }

    merged_column <- c(top_values, bot_values)
    merged_df[[col_name]] <- merged_column
  }

  return(merged_df)
}
