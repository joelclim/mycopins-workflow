library(dplyr)

get_lowest_evalue <- function(cluster_result_df, sequence_num) {
  cluster_seq_df <- cluster_result_df[
    cluster_result_df$sequence_num == sequence_num,]
  
  # no matches to evaluate
  if (nrow(cluster_seq_df) == 0) {
    return(c())
  }
  
  # keep the matches with the lowest evalue
  lowest_evalue <- min(cluster_seq_df$evalue)
  candidates_df <- cluster_seq_df[cluster_seq_df$evalue == lowest_evalue,]
  
  # remove the matches with a sciname that start with uncultured or fungal
  pattern <- "^uncultured|^fungal"   
  candidates_df <- candidates_df[
    !grepl(pattern, tolower(candidates_df$sciname)), ]
  
  scinames <- candidates_df$sciname
  
  # duplicates are allowed. occurrence number in search results matter.
  # More instance of the "species" means BLAST is certain. 
  # Otherwise, it is uncertain.
  return(scinames)
}


get_highest_percent_identity <- function(cluster_result_df, sequence_num) {
  cluster_seq_df <- cluster_result_df[
    cluster_result_df$sequence_num == sequence_num,]
  
  # no matches to evaluate
  if (nrow(cluster_seq_df) == 0) {
    return(c())
  }
  
  # keep the matches with the lowest evalue
  highest_percent_identity <- max(cluster_seq_df$percent_identity)
  candidates_df <- cluster_seq_df[
    cluster_seq_df$percent_identity == highest_percent_identity,]
  
  # remove the matches with a sciname that start with uncultured or fungal
  pattern <- "^uncultured|^fungal"   
  candidates_df <- candidates_df[!grepl(pattern, tolower(candidates_df$sciname)), ]
  
  scinames <- candidates_df$sciname
  
  # duplicates are allowed. occurrence number in search results matter.
  # More instance of the "species" means BLAST is certain. 
  # Otherwise, it is uncertain.
  return(scinames)
}


get_most_frequent <- function(cluster_result_df, all_candidates) {
  # Create a table to count occurrences of each value
  counts <- table(all_candidates)
  
  # Find the maximum count
  max_count <- max(counts)
  
  # Find the values with the maximum count
  most_frequent <- names(counts[counts == max_count])
  
  # Get the version of the most frequent with the smallest e-value.    
  subset_df <- subset(cluster_result_df, sciname %in% most_frequent)
  top_match <- subset_df[which.min(subset_df$evalue), ]
  
  return(top_match)
}


get_best_from_all_matches <- function(cluster_result_df) {
  # cluster_result_df was pre-sorted by increasing e-value and 
  # decreasing percent identity. 
  # Assume the first record is the top match.
  top_match <- cluster_result_df[1, ]
  
  # remove the matches with a sciname that start with uncultured or fungal
  pattern <- "^uncultured|^fungal"   
  candidates_df <- cluster_result_df[
    !grepl(pattern, tolower(cluster_result_df$sciname)), ]
  # Still remaining matches?
  if (nrow(candidates_df) > 1) {
    # Top of the remaining match has the same e-value and percent identity as
    # the first record of the top match.
    top_classified <- candidates_df[1, ]
    if (top_match$evalue == top_classified$evalue & 
        top_match$percent_identity == top_classified$percent_identity) {
      # Make the top classified the top match
      top_match <- top_classified
    }
  } else {
    # Otherwise, the first record is the top match
  }
  
  return(top_match)
}

get_best_match <- function(cluster_result_df) {
  match_by <- "LAST_RESORT"
  candidates_1 <- get_lowest_evalue(cluster_result_df, 1)
  candidates_2 <- get_lowest_evalue(cluster_result_df, 2)
  candidates_3 <- get_lowest_evalue(cluster_result_df, 3)
  all_candidates <- c(candidates_1, candidates_2, candidates_3)
  if (length(all_candidates) > 0) {
    best_match <- get_most_frequent(cluster_result_df, all_candidates)
    best_match$match_by <- "E-VALUE"
  } else {
    candidates_1 <- get_highest_percent_identity(cluster_result_df, 1)
    candidates_2 <- get_highest_percent_identity(cluster_result_df, 2)
    candidates_3 <- get_highest_percent_identity(cluster_result_df, 3)
    all_candidates <- c(candidates_1, candidates_2, candidates_3)
    if (length(all_candidates) > 0) {
      best_match <- get_most_frequent(cluster_result_df, all_candidates)
      best_match$match_by <- "PERCENT_IDENTITY"
    } else {
      best_match <- get_best_from_all_matches(cluster_result_df)
      best_match$match_by <- "BEST_EFFORT"
    }
  }
  
  return(best_match)
}

################################################################################
## 
## Reads the merged blast results and determine the top match for each cluster.
## A new CSV file is created containing the top match for each cluster.
##
################################################################################
top_blast_results_by_cluster <- function(consolidated_blast_df) {
  # Group together rows by cluster id. 
  # For each cluster group, sort based on cluster num, evalue exponent asc, 
  #   evalue coefficient desc, percent identity and hit num
  df <- consolidated_blast_df %>%
    group_by(cluster_id) %>%
    arrange(cluster_num, evalue_exponent, evalue_coefficient, 
            desc(percent_identity), hit_num)
  
  # Determine the cluster numbers and sort them in asc order
  cluster_nums <- sort(unique(df$cluster_num))
  
  # Create an empty data frame with the same structure as the merged data frame.
  top_match_df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  colnames(top_match_df) <- names(df)
  
  # For each cluster, determine the top match from the BLAST results.
  for (i in seq_along(cluster_nums)) {
    cluster_num <- cluster_nums[i]
    
    cluster_result_df <- df[df$cluster_num == cluster_num,]
    
    best_match <- get_best_match(cluster_result_df)
    
    print(paste0("Cluster ", cluster_num, " is ", best_match["sciname"], 
                 "; match by: ", best_match['match_by'],
                 "; e-value: ", best_match["evalue"], 
                 "; % identity: ", best_match["percent_identity"]))
    top_match_df <- rbind(top_match_df, best_match)
  }
  
  print("Sorting top match BLAST results by cluster number.")
  sorted_top_match_df <- top_match_df[order(top_match_df$cluster_num), ]
  
  print("Complete!")
  
  return(sorted_top_match_df)
}
