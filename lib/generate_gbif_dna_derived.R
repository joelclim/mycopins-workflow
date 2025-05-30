if (!require("readr")) install.packages("readr")
if (!require("tidyr")) install.packages("tidyr")

library(readr)
library(tidyr)

#
# Main function: generate_dna_derived
#
generate_dna_derived <- function(configuration) {

  get_samp_taxon_id <- function(source, taxid) {
    if (source == "BLAST") {
      return(paste0("NCBI:tx", taxid))
    }

    # UNITE
    tax_id <- strsplit(taxid, split="\\|")
    sh_id <- tax_id[[1]][1]
    seq_dataset <- tax_id[[1]][2]
    return(paste0("UNITE:", sh_id))
  }

  community_file <- configuration["mycopins_community_file"]
  mycopins.community <- read_csv(community_file, show_col_types = FALSE,
                                locale = locale(encoding = "UTF-8"))

  sample_cluster_file <- configuration["mycopins_sample_cluster_file"]
  mycopins.sample_cluster <- read_csv(sample_cluster_file, show_col_types = FALSE,
                                locale = locale(encoding = "UTF-8"))

  species <- unique(sort(names(mycopins.community)))

  gbif_dna_derived_df <- mycopins.sample_cluster %>%
    pivot_longer(
      cols = starts_with("Sequence"),
      names_to = "SeqRep",
      values_to = "DNA_sequence"
    ) %>%
    filter(
      (gbif_accepted_name %in% species) &
      !(gbif_accepted_name %in% c('Fungi', 'Fungus', 'mock')) &
      !(gbif.kingdom %in% c('Plantae')) &
      !is.na(DNA_sequence)
    ) %>%
    mutate(
      eventID = paste(Transect, Sample.Number, sep = "_"),
      occurrenceID = paste(paste(Transect, Sample.Number, sep = "_"), gbif.usage_key, sep=":"),
      samp_taxon_id = mapply(get_samp_taxon_id, source, taxid),
      lib_layout = configuration["libLayout"],
      target_gene = configuration["targetGene"],
      target_subfragment = configuration["targetSubfragment"],
      seq_meth = configuration["seqMethod"],
      otu_db = configuration["otuDB"],
      pcr_primer_name_forward = configuration["pcrPrimerNameFormward"],
      pcr_primer_name_reverse = configuration["pcrPrimerNameReverse"],
      pcr_primer_reference = configuration["pcrPrimerReference"]
    ) %>%
    select(eventID,
          occurrenceID,
          samp_name = Sample.Number,
          samp_taxon_id,
          lib_layout,
          target_gene,
          target_subfragment,
          seq_meth,
          otu_db,
          organism = gbif_accepted_name,
          pcr_primer_forward = Forward.Primer,
          pcr_primer_reverse = Reverse.Primer,
          pcr_primer_name_forward,
          pcr_primer_name_reverse,
          pcr_primer_reference,
          DNA_sequence
          )

  return(gbif_dna_derived_df)
}
