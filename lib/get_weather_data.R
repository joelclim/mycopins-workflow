if (!require("dplyr")) install.packages("dplyr")
if (!require("hms")) install.packages("hms")
if (!require("jsonlite")) install.packages("jsonlite")

library(dplyr)
library(hms)
library(jsonlite)

#
# Use to convert weather data in JSON format to a data frame.
# Set working directory where the weather data resides.
#
consolidate_weatherstack_results <- function(transect) {
  weather <- data.frame(
    transect = character(),
    date = character(),
    mintemp = numeric(),
    maxtemp = numeric(),
    totalsnow = numeric(),
    sunhour = numeric(),
    uv_index = numeric(),
    time = character(),
    temperature = numeric(),
    precip = numeric(),
    humidity = numeric(),
    feelslike = numeric()
  )

  json_files <- list.files(pattern = "\\.json$")
  for(json_file in json_files) {
    json_data <- fromJSON(json_file)
    historical <- json_data$historical
    for(day in historical) {
      hourly <- day$hourly
      for(index in rownames(hourly)) {
        row <- list(
          transect = transect,
          date = day$date,
          mintemp = day$mintemp,
          maxtemp = day$maxtemp,
          avgtemp = day$avgtemp,
          totalsnow = day$totalsnow,
          sunhour = day$sunhour,
          uv_index = day$uv_index,
          time = hourly[index, "time"],
          temperature = hourly[index, "temperature"],
          precip = hourly[index, "precip"],
          humidity = hourly[index, "humidity"],
          pressure = hourly[index, "pressure"], # atmospheric pressure
          heatindex = hourly[index, "heatindex"],
          dewpoint = hourly[index, "dewpoint"],
          feelslike = hourly[index, "feelslike"]
        )
        weather <- rbind(weather, as.data.frame(row, stringsAsFactors = FALSE))
      }
    }
    weather$date <- as.Date(weather$date)
  }

  return(weather)
}


get_weather_data <- function(weather_df, collection_date) {
  main_record <- weather_df %>%
    filter(date == collection_date) %>%
    slice(1) %>%
    select(mintemp, maxtemp, avgtemp, totalsnow, sunhour, uv_index)

  summary_record <- weather_df %>%
    filter(date == collection_date) %>%
    summarise(
      minprecip = min(precip),
      maxprecip = max(precip),
      avgprecip = mean(precip),
      minhumidity = min(humidity),
      maxhumidity = max(humidity),
      avghumidity = mean(humidity),
      minpressure = min(pressure),
      maxpressure = max(pressure),
      avgpressure = mean(pressure),
      minheatindex = min(heatindex),
      maxheatindex = max(heatindex),
      avgheatindex = mean(heatindex),
      mindewpoint = min(dewpoint),
      maxdewpoint = max(dewpoint),
      avgdewpoint = mean(dewpoint),
      minfeelslike = min(feelslike),
      maxfeelslike = max(feelslike),
      avgfeelslike = mean(feelslike)
    )

  return(cbind(main_record, summary_record))
}