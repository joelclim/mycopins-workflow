working_directory <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim"
data_directory <- paste0(working_directory, "/data/All")
environment_file <- paste0(data_directory, "/mycopins_environment_fungi.csv")
gbif_events_file <- paste0(data_directory, "/event.txt")

library(dplyr)
library(jsonlite)
library(readr)

mycopins.environment <- read_csv(environment_file, show_col_types = FALSE,
                         locale = locale(encoding = "UTF-8"))

location <- "Oulanka"
samplingProtocol <- "Shumskaya, M., Lorusso, N., Patel, U., Leigh, M., Somervuo, P., & Schigel, D. (2023). MycoPins: a metabarcoding-based method to monitor fungal colonization of fine woody debris. MycoKeys, 96, 77–95. https://doi.org/10.3897/mycokeys.96.101033"
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

gbif_events <- data.frame(
  eventID = character(),
  parentEventID = character(),
  samplingProtocol = character(),
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
  decimalLongitude = character()
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
      parentEventID <- paste0(transect, "_", date_id)
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
        eventID <- paste0(transect, "_", site)
        
        site_letter <- substr(site, nchar(site), nchar(site))
        
        dynamic_properties_df <- data.frame(woodType = get_wood_type(site_letter), 
                                            woodTexture = get_wood_texture(site_letter))
        weather_df <- data.frame(mintemp = env_record$mintemp,
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
                                 avgfeelslike = env_record$avgfeelslike)
        dynamic_properties_df$weather <- weather_df
        
        row <- list(
          eventID = eventID,
          parentEventID = parentEventID,
          samplingProtocol = samplingProtocol,
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
          dynamicProperties = gsub("\\[|\\]", "", toJSON(dynamic_properties_df))
        )
        gbif_events <- rbind(gbif_events, as.data.frame(row, stringsAsFactors = FALSE))
      }
    } else {
      print(paste("Transect", transect, "has nothing collected on", eventDate))
    }
  }
}

write.csv(gbif_events, gbif_events_file,
      row.names = FALSE, na = "", fileEncoding = "UTF-8")
