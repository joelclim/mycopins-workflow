library(readr)

scripts_path <- "C:/Users/Joel/work/kean-stme-2903-11/mycopins/scripts/"
source(paste0(scripts_path, "merge_add_counts.R"))

# Stage 2 - Taxonomy and Traits
source(paste0(scripts_path, "get_gbif_taxon.R"))
source(paste0(scripts_path, "get_fungal_traits.R"))

# Stage 3 - Fungi
source(paste0(scripts_path, "filter_for_fungi.R"))
source(paste0(scripts_path, "filter_for_wood_saprotroph.R"))

################################################################################
## Stage 1 - Merge (add)
################################################################################
working_folder_path <- "C:/Users/Joel/work/kean-stme-2903-11/mycopins/"

fungal_traits_path <- paste0(working_folder_path, "fungaltraits/")
fungal_traits_file <- paste0(fungal_traits_path, "fungal_traits_for_genera.csv")

fa23_path <- paste0(working_folder_path, "FA23/scata/")
fa23_counts_file <- paste0(fa23_path, "complete_counts.csv")
fa23_counts_df <- read.csv(fa23_counts_file, check.names=FALSE)
fa23_counts_df <- cbind(
  Batch = rep("FA23", nrow(fa23_counts_df)), fa23_counts_df
)

sp24_path <- paste0(working_folder_path, "SP24/scata/")
sp24_counts_file <- paste0(sp24_path, "complete_counts.csv")
sp24_counts_df <- read.csv(sp24_counts_file, check.names=FALSE)
sp24_counts_df <- cbind(
  Batch = rep("SP24", nrow(sp24_counts_df)), sp24_counts_df
)

merged_counts_df <- merge_add_counts(fa23_counts_df, sp24_counts_df)
merged_counts_file <- paste0(working_folder_path, 
                             "transectC_complete_counts.csv")
write.csv(merged_counts_df, file = merged_counts_file, row.names = FALSE)



################################################################################
## Stage 2 - Taxonomy and Traits
################################################################################
# Get gbif taxonomy for all identified organisms
sep_col_name <- "_Splitter_"
start_col <- which(names(merged_counts_df) == sep_col_name)
organisms <- names(merged_counts_df[,(start_col+1):ncol(merged_counts_df)])
gbif_taxon_df <- get_gbif_taxon(organisms)
gbif_taxon_file <- paste0(working_folder_path, 
                                   "transectC_gbif_taxon.csv")
write.csv(gbif_taxon_df, file = gbif_taxon_file, row.names = FALSE)

# Get fungal traits for genera of all identified organisms
fungal_traits_df <- read_csv(fungal_traits_file)
genus_traits_df <- get_fungal_traits(fungal_traits_df, gbif_taxon_df)
genus_traits_file <- paste0(working_folder_path, "transectC_genus_traits.csv")
write.csv(genus_traits_df, file = genus_traits_file, row.names = FALSE)


################################################################################
## Stage 3 - Fungi
################################################################################
# Complete counts containing only fungi
counts_fungi_df <- filter_for_fungi(merged_counts_df, gbif_taxon_df)
counts_fungi_file <- paste0(working_folder_path, "transectC_counts_fungi.csv")
write.csv(counts_fungi_df, file = counts_fungi_file, row.names = FALSE)

# Complete counts containing only wood saprotroph fungi
wood_saprotrophs_df <- get_wood_saprotrophs_traits(fungal_traits_df)
counts_wood_saprotroph_df <- filter_for_wood_saprotroph(
  merged_counts_df, gbif_taxon_df, wood_saprotrophs_df)
counts_wood_saprotroph_file <- paste0(working_folder_path, 
                                      "transectC_counts_wood_saprotroph.csv")
write.csv(counts_wood_saprotroph_df, 
          file = counts_wood_saprotroph_file, row.names = FALSE)
