if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")

library(dplyr)
library(readr)

clean_clusters <- function(in_all_clusters_file, out_all_clusters_file) {
  # Open the input file for reading
  in_con <- file(in_all_clusters_file, "r")

  # Open the output file for writing
  out_con <- file(out_all_clusters_file, "w")

  reference_count <- 0
  no_reference_count <- 0

  found_first_row <- FALSE
  # Read the input file line by line
  while (TRUE) {
    # Read a line from the input file
    line <- readLines(in_con, n = 1)

    # Check if the line is empty
    if (length(line) == 0) {
      break  # Exit the loop if end of file is reached
    }

    # Check if the line contains strings delimited by semicolon
    if (grepl(";", line)) {

      if (startsWith(line, "Cluster ID;Cluster Size;")) {
        # Header found
        found_first_row <- TRUE
      }
      if (found_first_row) {
        record <- line

        tokens <- strsplit(line, ";")[[1]]
        if (length(tokens) > 13) {
          reference_count <- reference_count + 1

          start_reference <- tokens[4]
          end_reference <- NULL
          for(i in seq_along(tokens)) {
            token <- tokens[i]
            if (i > 4 & !is.na(suppressWarnings(as.numeric(token)))) {
              end_reference <- token
              break
            }
          }
          left_end_index <- gregexpr(start_reference, line)[[1]][1]
          left_substring <- substring(line, 1, left_end_index-2)
          left_substring <- gsub(";", ",", left_substring)

          right_start_index <- gregexpr(end_reference, line)[[1]][1]
          right_substring <- substring(line, right_start_index, nchar(line))
          right_substring <- gsub(";", ",", right_substring)

          reference <- substring(line, left_end_index, right_start_index-1)
          reference <- gsub(",", "/", reference)

          record <- paste0(left_substring, ",", reference, ",", right_substring)
        } else {
          no_reference_count <- no_reference_count + 1
          record <- gsub(";", ",", record)
        }

        # Write the line to the output file
        cat(record, sep="\n", file = out_con)
      }
    }
  }

  # Close the input and output files
  close(in_con)
  close(out_con)

  suppressMessages(
    df <- read_csv(out_all_clusters_file,
                   show_col_types = FALSE,
                   locale = locale(encoding = "UTF-8"))
  )
  df <- select(df, -contains("..."))
  write.csv(df, file = out_all_clusters_file, row.names = FALSE, fileEncoding = "UTF-8")

  return(df)
}
