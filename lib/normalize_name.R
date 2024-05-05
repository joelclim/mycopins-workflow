normalize_name <- function(name) {
  new_name <- name

  # ignores
  if (new_name %in% c("uncultured organism")) {
    return(new_name)
  }

  # Remove " sp" or " sp." from the end of name, if exists.
  new_name <- gsub(" sp\\.?$", "", name)

  # Remove everything from " sp " and " sp. "
  # until the end of the name, if exists.
  new_name <- sub(" sp\\.? (.*)$", "", new_name)
  # "uncultured Trechispora" is "Teichospora"
  new_name <- gsub("^uncultured ", "", new_name)

  # uncultured fungus, fungal sp.
  if (new_name %in% c("fungus", "fungal")) {
    return("fungus")
  }
  if (startsWith(new_name, "fungal ")) {
    return ("fungus")
  }

  # mock 1, mock 2, mock 3
  if (startsWith(new_name, "mock ")) {
    return ("mock")
  }

  # TODO: foliar endophyte of Picea glauca?

  # Remove trailing white spaces
  return(trimws(new_name))
}
