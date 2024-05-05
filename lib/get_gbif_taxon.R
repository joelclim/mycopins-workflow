library("rgbif")

get_gbif_taxon <- function(organisms) {
  get_gbif_taxon <- function(organism) {
    print(paste("Retrieving gbif record for", organism))
    response <- name_backbone_verbose(name=organism, kingdom="fungi")
    taxon <- response$data
    if (taxon$matchType == "NONE") {
      return(c(
        matchType = taxon$matchType,
        scientificName = "",
        acceptedNameUsage = "",
        originalNameUsage = "",
        kingdom = "",
        phylum = "",
        class = "",
        order = "",
        family = "",
        genus = "",
        taxonKey = "",
        status = "",
        scientificNameAuthorship = ""
      ))
    }
    return(c(
      matchType = ifelse("matchType" %in% names(taxon), 
                         taxon$matchType, ""),
      scientificName = ifelse("scientificName" %in% names(taxon), 
                              taxon$scientificName, ""),
      acceptedNameUsage = ifelse("species" %in% names(taxon), 
                                 taxon$species, ""),
      originalNameUsage = ifelse("canonicalName" %in% names(taxon), 
                                 taxon$canonicalName, ""),
      kingdom = ifelse("kingdom" %in% names(taxon), taxon$kingdom, ""),
      phylum = ifelse("phylum" %in% names(taxon), taxon$phylum, ""),
      class = ifelse("class" %in% names(taxon), taxon$class, ""),
      order = ifelse("order" %in% names(taxon), taxon$order, ""),
      family = ifelse("family" %in% names(taxon), taxon$family, ""),
      genus = ifelse("genus" %in% names(taxon), taxon$genus, ""),
      taxonKey = ifelse("rank" %in% names(taxon), taxon$rank, ""),
      status = ifelse("status" %in% names(taxon), taxon$status, ""),
      scientificNameAuthorship = ifelse(
        "scientificName" %in% names(taxon) & 
        "canonicalName" %in% names(taxon), 
        trimws(substring(
          taxon$scientificName, 
          nchar(taxon$canonicalName)+1,
          nchar(taxon$scientificName))), "")
    ))
  }

  organism_gbif_df <- data.frame(organism = organisms)

  matchType <- character(nrow(organism_gbif_df))
  scientificName <- character(nrow(organism_gbif_df))
  acceptedNameUsage <- character(nrow(organism_gbif_df))
  originalNameUsage <- character(nrow(organism_gbif_df))
  kingdom <- character(nrow(organism_gbif_df))
  phylum <- character(nrow(organism_gbif_df))
  class <- character(nrow(organism_gbif_df))
  order <- character(nrow(organism_gbif_df))
  family <- character(nrow(organism_gbif_df))
  genus <- character(nrow(organism_gbif_df))
  taxonKey <- character(nrow(organism_gbif_df))
  status <- character(nrow(organism_gbif_df))
  scientificNameAuthorship <- character(nrow(organism_gbif_df))

  for (i in 1:nrow(organism_gbif_df)) {
    taxon <- get_gbif_taxon(organism_gbif_df[i, ])

    matchType[i] <- taxon['matchType']
    scientificName[i] <- taxon['scientificName']
    acceptedNameUsage[i] <- taxon['acceptedNameUsage']
    originalNameUsage[i] <- taxon['originalNameUsage']
    kingdom[i] <- taxon['kingdom']
    phylum[i] <- taxon['phylum']
    class[i] <- taxon['class']
    order[i] <- taxon['order']
    family[i] <- taxon['family']
    genus[i] <- taxon['genus']
    taxonKey[i] <- taxon['taxonKey']
    status[i] <- taxon['status']
    scientificNameAuthorship[i] <- taxon['scientificNameAuthorship']
  }
  
  organism_gbif_df$matchType <- matchType
  organism_gbif_df$scientificName <- scientificName
  organism_gbif_df$acceptedNameUsage <- acceptedNameUsage
  organism_gbif_df$originalNameUsage <- originalNameUsage
  organism_gbif_df$kingdom <- kingdom
  organism_gbif_df$phylum <- phylum
  organism_gbif_df$class <- class
  organism_gbif_df$order <- order
  organism_gbif_df$family <- family
  organism_gbif_df$genus <- genus
  organism_gbif_df$taxonKey <- taxonKey
  organism_gbif_df$status <- status
  organism_gbif_df$scientificNameAuthorship <- scientificNameAuthorship
  
  return(organism_gbif_df)
}
