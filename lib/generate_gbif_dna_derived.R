library(dplyr)
library(readr)

#
# Main function: generate_dna_derived
#
generate_dna_derived <- function(configuration) {
  tag_sequence_file <- configuration["mycopins_tag_sequence_file"]
  mycopins.tag_sequence <- read_csv(tag_sequence_file, show_col_types = FALSE,
                                locale = locale(encoding = "UTF-8"))

  gbif_dna_derived_df <- mycopins.tag_sequence %>%
    pivot_longer(
      cols = starts_with("Sequence"),
      names_to = "SeqRep",
      values_to = "DNA_sequence"
    ) %>%
    filter(
      !(gbif_accepted_name %in% c("Fungus", "fungi")) &
      !is.na(DNA_sequence)
    ) %>%
    mutate(
      eventID = paste(Transect, Sample.Number, sep = "_"),
      occurrenceID = paste(paste(Transect, Sample.Number, sep = "_"), gbif.usage_key, sep=":"),
      samp_taxon_id = paste0("https://www.gbif.org/species/", gbif.usage_key)
    ) %>%
    select(eventID,
           occurrenceID,
           samp_name = Sample.Number,
           samp_taxon_id,
           organism = gbif_accepted_name,
           pcr_primer_forward = Forward.Primer,
           pcr_primer_reverse = Reverse.Primer,
           DNA_sequence
          )

  return(gbif_dna_derived_df)
}
