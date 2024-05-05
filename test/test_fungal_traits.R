install.packages("devtools")
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("traitecoevo/fungaltraits")

library(fungaltraits)
df <- fungal_traits()

working_folder_path <- "C:/Users/Joel/work/kean-stme-2903-11/mycopins/"
mycopins_folder_path <- paste0(working_folder_path, "SP24/scata/")
fungal_traits_path <- paste0(working_folder_path, "fungaltraits/")

library(readr)
fungal_traits_file <- paste0(fungal_traits_path, "fungal_traits_for_genera.csv")
fungal_traits_df <- read_csv(fungal_traits_file)
length(unique(fungal_traits_df$GENUS))

organism_gbif_taxon_file <- paste0(mycopins_folder_path, "organism_gbif_taxon.csv")
mycopins_df <- read_csv(organism_gbif_taxon_file)
has_genus_df <- mycopins_df[!is.na(mycopins_df$genus), ]
mycopin_genus <- unique(sort(has_genus_df$genus))
mycopin_traits_df <- fungal_traits_df[fungal_traits_df$GENUS %in% mycopin_genus,]
