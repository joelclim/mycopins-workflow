if (!require("readr")) install.packages("readr")

library(readr)

batch_cluster_sequences <- function(cleaned_clusters_df,
                                    blast_folder_path,
                                    batch_file_prefix = "blast_cluster_sequences",
                                    max_batch_size = 100) {
  batch_path_prefix <- paste0(blast_folder_path, "blast_cluster_sequences")

  no_reference_df <- cleaned_clusters_df[is.na(cleaned_clusters_df$Reference), ]
  search_df <- no_reference_df

  lines <- NULL
  batch <- 0
  count <- 0
  for (i in seq_len(nrow(search_df))) {
    row <- search_df[i, ]

    lines <- c(lines, paste0("> ", row["Cluster ID"], "_1"))
    lines <- c(lines, row$Sequence1)
    if (!is.na(row$Sequence2)) {
      lines <- c(lines, paste0("> ", row["Cluster ID"], "_2"))
      lines <- c(lines, row$Sequence2)
    }
    if (!is.na(row$Sequence3)) {
      lines <- c(lines, paste0("> ", row["Cluster ID"], "_3"))
      lines <- c(lines, row$Sequence3)
    }

    count <- count + 1

    if (max_batch_size == count) {
      batch <- batch + 1

      # write a batch of sequences to a file.
      out_file <- paste0(batch_path_prefix, batch, ".fasta")
      file_con <- file(out_file)
      writeLines(lines, file_con)
      close(file_con)

      lines <- NULL
      count <- 0
    }
  }

  # write the last remaining batch of sequences to a file.
  batch <- batch + 1
  out_file <- paste0(batch_path_prefix, batch, ".fasta")
  file_con <- file(out_file)
  writeLines(lines, file_con)
  close(file_con)
}