get_fungal_traits <- function(fungal_traits_df, organism_df) {
  has_genus_df <- organism_df[!is.na(organism_df$gbif.genus), ]
  genus <- unique(sort(has_genus_df$gbif.genus))

  genus_traits_df <- fungal_traits_df[fungal_traits_df$GENUS %in% genus,]
  return(genus_traits_df)
}

get_wood_saprotrophs_traits <- function(traits_df) {
  primary_traits_df <- traits_df[traits_df$primary_lifestyle
                                 == "wood_saprotroph" &
                                   !is.na(traits_df$primary_lifestyle), ]
  wood_saprotroph_genera <- primary_traits_df$GENUS
  traits_df[traits_df$Secondary_lifestyle == "wood_saprotroph", ]
  secondary_traits_df <- traits_df[traits_df$Secondary_lifestyle
                                   == "wood_saprotroph" &
                                     !is.na(traits_df$Secondary_lifestyle), ]
  wood_saprotroph_genera <- sort(unique(c(wood_saprotroph_genera,
                                          secondary_traits_df$GENUS)))

  wood_saprotrophs_traits_df <- rbind(primary_traits_df,
                                      secondary_traits_df)

  return(wood_saprotrophs_traits_df)
}