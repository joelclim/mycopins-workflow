working_directory <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim/mycopins"
setwd(working_directory)

batch_name <- "FA23-transectC"
scata_dataset_name <- "FA23-C-Mycopins"
scata_job_name <- "scata6398"

#batch_name <- "SP24-transectC-1"
#scata_dataset_name <- "SP24-RFI-TranC \\(Full\\)"
#scata_job_name <- "scata6397"

#batch_name <- "SP24-transectC-2"
#scata_dataset_name <- "SP24-transectC-2"
#scata_job_name <- "scata6466"

################################################################################
# Data sources and targets (convention-over-configuration approach)
################################################################################
workflow_path <- "./workflow/"
lib_path <- paste0(workflow_path, "lib/")

reference_data_directory <- "./reference-data/"
fungal_traits_directory <- paste0(reference_data_directory, "fungaltraits/")
fungal_traits_file <- paste0(fungal_traits_directory, "fungal_traits_for_genera.csv")

data_directory <- "./data/"
batch_directory <- paste0(data_directory, batch_name, "/")
scata_directory <- paste0(batch_directory, "scata/")
scata_job_directory <- paste0(scata_directory, scata_job_name, "/")
scata_counts_file <- paste0(scata_job_directory, "all_tag_by_cluster_counts.txt")
scata_clusters_file <- paste0(scata_job_directory, "all_clusters_", scata_job_name, ".txt")

blast_directory <- paste0(batch_directory, "blast/")

features_file_name <- paste0(batch_name, "-features.csv")
features_file <- paste0(batch_directory, features_file_name)
################################################################################
## IMPORTS
################################################################################
library(readr)

# Stage 1 - Clean
source(paste0(lib_path, "clean_counts.R"))
source(paste0(lib_path, "clean_clusters.R"))
source(paste0(lib_path, "batch_cluster_sequences.R"))

# Stage 2 - Identify
source(paste0(lib_path, "consolidate_blast_results.R"))
source(paste0(lib_path, "top_blast_results_by_cluster.R"))
source(paste0(lib_path, "normalize_name.R"))
source(paste0(lib_path, "get_gbif_accepted_name.R"))
source(paste0(lib_path, "apply_top_match_to_clusters.R"))
source(paste0(lib_path, "apply_match_to_counts.R"))

# Stage 3 - Taxonomy and Traits
source(paste0(lib_path, "get_gbif_taxon.R"))
source(paste0(lib_path, "get_fungal_traits.R"))

# Stage 4 - Fungi
source(paste0(lib_path, "filter_for_fungi.R"))
source(paste0(lib_path, "filter_for_wood_saprotroph.R"))


################################################################################
## WORKFLOW
################################################################################

################################################################################
## Stage 1 - Clean
################################################################################

# Clean counts
cleaned_counts_df <- clean_counts(scata_counts_file, features_file, scata_dataset_name)
cleaned_counts_file <- paste0(batch_directory, "cleaned_counts.csv")
write.csv(cleaned_counts_df, file = cleaned_counts_file, row.names = FALSE)

# Clean clusters
cleaned_clusters_file <- paste0(batch_directory, "cleaned_clusters.csv")
cleaned_clusters_df <- clean_clusters(scata_clusters_file, cleaned_clusters_file)

# Blast cluster sequences
batch_cluster_sequences(cleaned_clusters_df, blast_directory)

################################################################################
# 1. Go to NCBI (blastn): https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn
# 2. Choose File > (batch clusters sequences#.fasta)
# 3. Algorithm Parameters > General Parameters > Max target sequences = 10
# 4. Press the BLAST button.
# 5. Download All > Single-file JSON > Save > (blast_directory)/(results).json
################################################################################


################################################################################
## Stage 2 - Identify
################################################################################

# Consolidate BLAST results
consolidated_blast_df <- consolidate_blast_results(blast_directory, scata_job_name)
consolidated_blast_file <- paste0(batch_directory, "blast_consolidated.csv")
write.csv(consolidated_blast_df, file = consolidated_blast_file, row.names = FALSE)

# Top BLAST results by cluster
#consolidated_blast_df <- read_csv(consolidated_blast_file)
top_match_df <- top_blast_results_by_cluster(consolidated_blast_df)
top_match_file <- paste0(batch_directory, "blast_top_match.csv")
write.csv(top_match_df, file = top_match_file, row.names = FALSE)

# Apply top match to clusters
#top_match_df <- read_csv(top_match_file)
cleaned_clusters_df <- read_csv(cleaned_clusters_file)
complete_clusters_df <- apply_top_match_to_clusters(top_match_df, cleaned_clusters_df)
complete_clusters_file <- paste0(batch_directory, "complete_clusters.csv")
write.csv(complete_clusters_df, file = complete_clusters_file, row.names = FALSE)

# Apply match to counts
#cleaned_counts_df <- read_csv(cleaned_counts_file)
#complete_clusters_df <- read_csv(complete_clusters_file)
complete_counts_df <- apply_match_to_counts(cleaned_counts_df, complete_clusters_df)
complete_counts_file <- paste0(batch_directory, "complete_counts.csv")
write.csv(complete_counts_df, file = complete_counts_file, row.names = FALSE)


################################################################################
## Stage 3 - Taxonomy and Traits
################################################################################
# Get gbif taxonomy for all identified organisms
#complete_clusters_df <- read_csv(complete_clusters_file)
organisms <- unique(sort(complete_clusters_df$normalized_name))
gbif_taxon_df <- get_gbif_taxon(organisms)
gbif_taxon_file <- paste0(batch_directory, "batch_gbif_taxon.csv")
write.csv(gbif_taxon_df, file = gbif_taxon_file, row.names = FALSE)

# Get fungal traits for genera of all identified organisms
#gbif_taxon_df <- read_csv(batch_gbif_taxon_file)
fungal_traits_df <- read_csv(fungal_traits_file)
genus_traits_df <- get_fungal_traits(fungal_traits_df, gbif_taxon_df)
genus_traits_file <- paste0(batch_directory, "batch_genus_traits.csv")
write.csv(genus_traits_df, file = genus_traits_file, row.names = FALSE)


################################################################################
## Stage 4 - Fungi
################################################################################
# Complete counts containing only fungi
counts_fungi_df <- filter_for_fungi(complete_counts_df, gbif_taxon_df)
counts_fungi_file <- paste0(batch_directory, "complete_counts_fungi.csv")
write.csv(counts_fungi_df, file = counts_fungi_file, row.names = FALSE)

# Complete counts containing only wood saprotroph fungi
wood_saprotrophs_df <- get_wood_saprotrophs_traits(fungal_traits_df)
counts_wood_saprotroph_df <- filter_for_wood_saprotroph(
  complete_counts_df, gbif_taxon_df, wood_saprotrophs_df)
counts_wood_saprotroph_file <- paste0(batch_directory, 
                                      "complete_counts_wood_saprotroph.csv")
write.csv(counts_wood_saprotroph_df, 
          file = counts_wood_saprotroph_file, row.names = FALSE)
