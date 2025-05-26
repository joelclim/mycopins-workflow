library(dplyr)
library(readr)

# Load the dataset
occurrence_df <- read_csv("occurrence-all-fungi.txt")

occurrence_df <- occurrence_df %>%
  filter(occurrenceStatus == "present") %>%
  select(originalNameUsage, organismQuantity, taxonRank) %>%
  arrange(originalNameUsage) 

# Count occurrences per originalNameUsage and taxonRank
occurrence_df <- occurrence_df %>%
  group_by(originalNameUsage, taxonRank) %>%
  summarise(occurrenceCount = n(), .groups = "drop")

# Get the count and sum of total_quantity per taxonRank
occurrence_df <- occurrence_df %>%
  group_by(taxonRank) %>%
  summarise(totalNumberOfTaxa = n(), numberOfOTUsIdentifiedAtThisLevel = sum(occurrenceCount), .groups = "drop")

# Calculate percentage of sum_counts relative to overall numberOfOTUsIdentifiedAtThisLevel
total_sum_counts <- sum(occurrence_df$numberOfOTUsIdentifiedAtThisLevel)
occurrence_df <- occurrence_df %>%
  mutate(percentageOfOTUsIdentifiedAtThisLevel = (numberOfOTUsIdentifiedAtThisLevel / total_sum_counts) * 100)

# Sort based on taxonRank order
rank_order <- c("KINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")
occurrence_df$taxonRank <- factor(occurrence_df$taxonRank, levels = rank_order)
occurrence_df <- occurrence_df %>% arrange(taxonRank)

write_csv(occurrence_df, "taxonomic_coverage.csv")
