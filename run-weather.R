workflow_directory <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim/mycopins-workflow"
source(paste0(workflow_directory, "/lib/get_weather_data.R"))

weather_data_directory <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim/weather-data"

transect <- "AandB"
working_directory <- paste0(weather_data_directory, "/transect", transect)
setwd(working_directory)

weather <- consolidate_weatherstack_results("A")
write_csv(weather, paste0(output_directory, "/weather-A.csv"))
weather <- consolidate_weatherstack_results("B")
write_csv(weather, paste0(output_directory, "/weather-B.csv"))

transect <- "C"
working_directory <- paste0(weather_data_directory, "/transect", transect)
setwd(working_directory)
weather <- consolidate_weatherstack_results(transect)
write_csv(weather, paste0(working_directory, "/weather.csv"))

transect <- "D"
working_directory <- paste0(weather_data_directory, "/transect", transect)
setwd(working_directory)
weather <- consolidate_weatherstack_results(transect)
write_csv(weather, paste0(working_directory, "/weather.csv"))
