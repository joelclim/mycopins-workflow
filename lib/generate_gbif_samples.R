working_directory <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim"
data_directory <- paste0(working_directory, "/data/All")
environment_file <- paste0(data_directory, "/mycopins_environment_fungi.csv")
gbif_samples_file <- paste0(data_directory, "/Samples.tsv")

library(dplyr)
library(jsonlite)
library(readr)

mycopins.environment <- read_csv(environment_file, show_col_types = FALSE,
                                 locale = locale(encoding = "UTF-8"))

location <- "Oulanka"
locationID <- "https://www.geonames.org/12226273"
countryCode <- "FI"
country <- "Finland"
stateProvince <- "North Ostrobothnia"
municipality <- "Kuusamo"

transects <- c("A", "C")
habitats <- c(
  "A" = "boreal forest area protected from grazing by reindeers; wooden pins buried in soil.",
  "C" = "mixed broadleaf forest accessed by random visitors; wooden pins buried in soil."
)
latitude <- c(
  "A" = 66.367,
  "C" = 66.376
)
longitude <- c(
  "A" = 29.533,
  "C" = 29.313
)

get_wood_type <- function(site) {
  if (site %in% c("A", "B")) {
    return("Pine")
  }
  if (site %in% c("C", "D")) {
    return("Birch")
  }
  # site %in% c("E", "F")
  return ("Spruce")
}

get_wood_texture <- function(site) {
  if (site %in% c("A", "B", "E", "F")) {
    return("Softwood")
  }
  # site %in% c("C", "D")
  return ("Hardwood")
}

gbif_samples <- data.frame(
  id = character(),
  samplingEffort = character(),
  eventDate = character(),
  habitat = character(),
  fieldNotes = character(),
  locationID = character(),
  countryCode = character(),
  country = character(),
  stateProvince = character(),
  municipality = character(),
  decimalLatitude = character(),
  decimalLongitude = character(),
  woodType = character(),
  woodTexture = character(),
  mintemp = numeric(),
  maxtemp = numeric(),
  avgtemp = numeric(),
  minhumidity = numeric(),
  maxhumidity = numeric(),
  avghumidity = numeric(),
  minpressure = numeric(),
  maxpressure = numeric(),
  avgpressure = numeric(),
  minheatindex = numeric(),
  maxheatindex = numeric(),
  avgheatindex = numeric(),
  mindewpoint = numeric(),
  maxdewpoint = numeric(),
  avgdewpoint = numeric(),
  minfeelslike = numeric(),
  maxfeelslike = numeric(),
  avgfeelslike = numeric()
)

for (i in 1:length(transects)) {
  transect <- transects[i]
  collection_dates <- mycopins.environment$Date.Collected[
    mycopins.environment$Transect == transect
  ]
  collection_dates <- as.Date(collection_dates, "%m/%d/%Y")
  collection_dates <- sort(unique(collection_dates))
  for (j in 1:length(collection_dates)) {
    eventDate <- collection_dates[j]
    print(paste("Transect:", transect, "; Date.Collected:", eventDate))
    
    date_id <- format(eventDate, "%Y_%b_%d")
    env_record <- mycopins.environment %>%
      filter(as.Date(Date.Collected, "%m/%d/%Y") == eventDate 
             & Transect == transect) %>%
      slice(1)
    if (nrow(env_record) == 1) { # exists
      samplingEffort <- paste(env_record$Days.Elapsed, "days")
      fieldNotes <- env_record$Season
      habitat <- habitats[transect]
      decimalLatitude <- latitude[transect]
      decimalLongitude <- longitude[transect]
      
      sites <- mycopins.environment %>%
        filter(as.Date(Date.Collected, "%m/%d/%Y") == eventDate 
               & Transect == transect) %>%
        arrange(Sample.Number) %>%
        pull(Sample.Number)
      
      for (k in 1:length(sites)) {
        site <- sites[k]
        site_letter <- substr(site, nchar(site), nchar(site))

        row <- list(
          id = site,
          samplingEffort = samplingEffort,
          eventDate = eventDate,
          habitat = habitat,
          fieldNotes = fieldNotes,
          locationID = locationID,
          countryCode = countryCode,
          country = country,
          stateProvince = stateProvince,
          municipality = municipality,
          decimalLatitude = decimalLatitude,
          decimalLongitude = decimalLongitude,
          woodType = get_wood_type(site_letter),
          woodTexture = get_wood_texture(site_letter),
          mintemp = env_record$mintemp,
          maxtemp = env_record$maxtemp,
          avgtemp = env_record$avgtemp,
          minhumidity = env_record$minhumidity,
          maxhumidity = env_record$maxhumidity,
          avghumidity = env_record$avghumidity,
          minpressure = env_record$minpressure,
          maxpressure = env_record$maxpressure,
          avgpressure = env_record$avgpressure,
          minheatindex = env_record$minheatindex,
          maxheatindex = env_record$maxheatindex,
          avgheatindex = env_record$avgheatindex,
          mindewpoint = env_record$mindewpoint,
          maxdewpoint = env_record$maxdewpoint,
          avgdewpoint = env_record$avgdewpoint,
          minfeelslike = env_record$minfeelslike,
          maxfeelslike = env_record$maxfeelslike,
          avgfeelslike = env_record$avgfeelslike
        )
        gbif_samples <- rbind(gbif_samples, as.data.frame(row, stringsAsFactors = FALSE))
      }
    } else {
      print(paste("Transect", transect, "has nothing collected on", eventDate))
    }
  }
}

write.table(gbif_samples, gbif_samples_file, sep = "\t", 
            row.names = FALSE, na = "", fileEncoding = "UTF-8")
