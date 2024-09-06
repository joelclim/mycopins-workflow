library(dplyr)
library(jsonlite)
library(readr)

#
# Requires get_transects(), get_habitat(transect), get_latitude(transect), get_longitude(transect)
#

create_gbif_event_dynamic_properties <- function(site_letter, organism, env_record) {
  dynamic_properties_df <- data.frame(woodType = get_wood_type(site_letter),
                                      woodTexture = get_wood_texture(site_letter))
  weather_df <- data.frame(
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
  dynamic_properties_df$weather <- weather_df

  return(dynamic_properties_df)
}

generate_gbif_event <- function(configuration) {
  environment_file <- configuration["mycopins_environment_file"]

  location <- configuration["location"]
  samplingProtocol <- configuration["samplingProtocol"]
  locationID <- configuration["locationID"]
  countryCode <- configuration["countryCode"]
  country <- configuration["country"]
  stateProvince <- configuration["stateProvince"]
  municipality <- configuration["municipality"]

  transects <- get_transects()

  mycopins.environment <- read_csv(environment_file, show_col_types = FALSE,
                                   locale = locale(encoding = "UTF-8"))

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
        decimalLatitude <- get_latitude(transect)
        decimalLongitude <- get_longitude(transect)

        sites <- mycopins.environment %>%
          filter(as.Date(Date.Collected, "%m/%d/%Y") == eventDate
                & Transect == transect) %>%
          arrange(Sample.Number) %>%
          pull(Sample.Number)

        for (k in 1:length(sites)) {
          site <- sites[k]
          eventID <- paste0(transect, "_", site)

          site_letter <- substr(site, nchar(site), nchar(site))
          wood_type <- get_wood_type(site_letter)
          wood_texture <- get_wood_texture(site_letter)
          habitat <- paste0(wood_type, "(", wood_texture, ") dowel ", get_habitat(transect))

          dynamic_properties_df <- create_gbif_event_dynamic_properties(site_letter, organism, env_record)

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

  return(gbif_events)
}
