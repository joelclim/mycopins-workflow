if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
library(dplyr)
library(tidyr)

create_tag_sequence_dataset <- function(cleaned_counts_df, complete_clusters_df, mycopins_organisms_df) {
  # Reshape data to a long format for scata cluster ids
  long_cleaned_counts_df <- cleaned_counts_df %>%
    pivot_longer(
      cols = starts_with("scata"),
      names_to = "Cluster.ID",
      values_to = "Count"
    ) %>%
    filter(Count > 0)

  tag_sequence_df <- long_cleaned_counts_df %>%
    inner_join(complete_clusters_df, by = "Cluster.ID") %>%
    select(Transect, Sample.Number,     # clean_counts_df
      Cluster.ID, Count,
      Forward.Primer, Reverse.Primer,
      Sequence1, Sequence2, Sequence3,
      taxid, gbif_accepted_name         # complete_clusters_df
    )

  tag_sequence_with_gbif_taxon_id_df <- tag_sequence_df %>%
    inner_join(mycopins_organisms_df, by = c("gbif_accepted_name"="organism")) %>%
    select(Transect, Sample.Number,                   # clean_counts_df
      Cluster.ID, Count,
      Forward.Primer, Reverse.Primer,
      Sequence1, Sequence2, Sequence3,
      taxid, gbif_accepted_name,                      # complete_clusters_df
      match.source, match.tax_id, gbif.usage_key      # mycopins_organisms_df
    )

  return(tag_sequence_with_gbif_taxon_id_df %>%
    filter(!is.na(match.source)))
}
