
library(ncdf4)
library(lubridate)
library(circular)
library(dplyr)
library(tidyverse)
library(tidyr)

# This script shows the necessary steps to obtain the data for the San Juan de Fuca
# analysis. There are two main steps: (1) downloading the raw data;
# (2) filtering the raw data to obtain data only for zones A, B, C, and D with
# the weighted average daily mean

## Obtaining raw data files for the full west coast

# ID for west coast
data_ID <- "HFR/USWC/2km/hourly/RTV/HFRADAR_US_West_Coast_2km_Resolution_Hourly_RTV_best.ncd"

# Data url
data_url <- paste0("http://hfrnet-tds.ucsd.edu/thredds/dodsC/", data_ID)

# Read data structure
data <- nc_open(data_url)

# Retrieve latitudes, longitudes, and time
lat <- ncvar_get(data, "lat")
lon <- ncvar_get(data, "lon")
tim <- date("2012-01-01 00:00:00 UTC") + hours(ncvar_get(data, "time"))

# Trim download region -- server has a limit of 500Mb request
# Check (lat, lon) coordinates at https://cordc.ucsd.edu/projects/mapping/maps/
# to focus on a given region
loc <- "Strait_of_Juan_de_Fuca"
begin_lat <- 48.1
end_lat <- 48.5
begin_lon <- -124.5
end_lon <- -123

# Get data function -- takes info from global environment and writes on it
get_data <- function() {

  # Obtain the data indexes associated to the given information
  begin_lat_ind <<- which(lat >= begin_lat)[1]
  begin_lon_ind <<- which(lon >= begin_lon)[1]
  begin_tim_ind <<- which(tim >= begin_tim)[1]
  end_lat_ind <<- which(lat <= end_lat)
  end_lon_ind <<- which(lon <= end_lon)
  end_tim_ind <<- which(tim <= end_tim)
  end_lat_ind <<- end_lat_ind[length(end_lat_ind)]
  end_lon_ind <<- end_lon_ind[length(end_lon_ind)]
  end_tim_ind <<- end_tim_ind[length(end_tim_ind)]

  # Download sizes
  l_lat <<- end_lat_ind - begin_lat_ind + 1
  l_lon <<- end_lon_ind - begin_lon_ind + 1
  l_tim <<- end_tim_ind - begin_tim_ind + 1
  stopifnot(l_lat > 0)
  stopifnot(l_lon > 0)
  stopifnot(l_tim > 0)

  # Upper bound on the size of one out of two objects to be downloaded
  cat("Data size:", format(object.size(rnorm(l_lat * l_lon * l_tim)),
                           units = "Mb"))

  # Download (u, v), being:
  # u (m/s) = surface_eastward_sea_water_velocity
  # v (m/s) = surface_northward_sea_water_velocity
  u <<- ncvar_get(data, "u",
                  start = c(begin_lon_ind, begin_lat_ind, begin_tim_ind),
                  count = c(l_lon, l_lat, l_tim))
  v <<- ncvar_get(data, "v",
                  start = c(begin_lon_ind, begin_lat_ind, begin_tim_ind),
                  count = c(l_lon, l_lat, l_tim))

}

# Download data in a monthly loop to avoid surpassing 500Mb limit
for (year in 2021:2021) {
  
  month <- 1
  if (year == 2020) {
    month <- 6
  }
  while ( ( year < 2022 && month <= 12 ) || (year == 2022 && month <= 6) ) {

    # Show progress
    cat("\n", year, "-", month, "\n")

    begin_tim <- date(paste(toString(year), toString(month),
                            "01 00:00:00 UTC", sep = "-"))
    if (month != 12) {

      end_tim <- date(paste(toString(year), toString(month + 1),
                            "01 00:00:00 UTC", sep = "-"))

    } else {

      end_tim <- date(paste(toString(year + 1), toString(1),
                            "01 00:00:00 UTC", sep = "-"))

    }

    # Download data
    get_data()

    # Get data dimensions
    lat_length <- dim(u)[1]
    lon_length <- dim(u)[2]
    time_length <- dim(u)[3]

    # Find the indexes associated to dimensions
    lat_aux <- begin_lat_ind:(begin_lat_ind + lat_length - 1)
    lon_aux <- begin_lon_ind:(begin_lon_ind + lon_length - 1)
    time_aux <- begin_tim_ind:(begin_tim_ind + time_length - 1)

    # Join all the cases
    join <- merge(x = lat[lat_aux], y = lon[lon_aux])
    final_dataframe <- merge(x = join, y = tim[time_aux], by = NULL)
    colnames(final_dataframe) <- c("lat", "lon", "time")

    # Add the velocities and save the data
    final_dataframe$u <- c(u)
    final_dataframe$v <- c(v)
    save(final_dataframe, file = paste(loc, "_", toString(year), "_",
                                       toString(month), ".RData", sep = ""))
    
    month <- month + 1
  }

}

## Obtaining daily data for zones A, B, C, and D

# Function that reads over all the files in the directory containing the raw
# data and return the records inside a given area delimited by longitude
# and latitude
extract_data <- function(begin_lat, end_lat, begin_lon, end_lon) {
  
  # Retrieve monthly data
  monthly_data <- lapply(files, function(x) {
    
    # Load raw data file
    load(x)
    results <- filter(final_dataframe, lat > begin_lat, lat < end_lat,
                      lon > begin_lon, lon < end_lon)
    
    # Calculate directions and speed
    results$d <- atan2(x = results$u, y = results$v)
    results$speed <- sqrt(results$u^2 + results$v^2)
    return(results)
    
  })
  
  # Merge available data
  total_data <- data.frame()
  for (i in seq_along(monthly_data)) {
    
    total_data <- rbind(total_data, monthly_data[[i]])
    
  }
  return(total_data)
  
}

# Function that calculates the speed-weighted average of the directions over
# periods of certain number of hours of duration
extract_datetime_theta <- function(results, hours) {
  # Order data by ascending time
  results <- results[order(as.Date(results$time)),]
  # Remove NA's
  results <- results[complete.cases(results$u, results$v, results$d, results$speed),]
  results %>% 
    # Assigning an interval of three hours to each observation
    mutate(
      interval = as.POSIXct(cut(time, breaks = paste(hours, "hours")))
    ) %>% 
    # Grouping by the variable interval
    group_by(interval) %>%
    mutate(
      weights = speed / sum(speed),
      theta=circular:::WeightedMeanCircularRad(w = weights, x = d)
    ) %>% 
    ungroup %>% 
    select(interval, theta) %>% 
    unique %>%
    data.frame
}


# Zone A
begin_lat_A <- 48.1
end_lat_A <- 48.5
begin_lon_A <- -124.5
end_lon_A <- -124.3125

# Zone B
begin_lat_B <- 48.1
end_lat_B <- 48.5
begin_lon_B <- -124.3125
end_lon_B <- -124.25

# Zone C
begin_lat_C <- 48.1
end_lat_C <- 48.5
begin_lon_C <- -124.25
end_lon_C <- -124.1875

# Zone D
begin_lat_D <- 48.1
end_lat_D <- 48.5
begin_lon_D <- -124.1875
end_lon_D <- -124.125

# List individual RDatas
files <- list.files(pattern = "*.RData", full.names = TRUE, recursive = FALSE)

# Obtain the results for every area and their corresponding daily directions
results_A <- extract_data(begin_lat_A, end_lat_A, begin_lon_A, end_lon_A)
A <- extract_datetime_theta(results_A, hours = 3)

results_B <- extract_data(begin_lat_B, end_lat_B, begin_lon_B, end_lon_B)
B <- extract_datetime_theta(results_B, hours = 3)

results_C <- extract_data(begin_lat_C, end_lat_C, begin_lon_C, end_lon_C)
C <- extract_datetime_theta(results_C, hours = 3)

results_D <- extract_data(begin_lat_D, end_lat_D, begin_lon_D, end_lon_D)
D <- extract_datetime_theta(results_D, hours = 3)

df_list <- list(A, B, C, D)

data <- df_list %>% 
  reduce(full_join, by='interval') %>% 
  rename(zone_A = theta.x,
         zone_B = theta.y,
         zone_C = theta.x.x,
         zone_D = theta.y.y)
data <- data[complete.cases(data$zone_A, data$zone_B, data$zone_C, data$zone_D),]

sanjuanfuca <- array(NA, dim = c(nrow(data), 2, 4))
sanjuanfuca[ , , 1] = DirStats::to_cir(data$zone_A)
sanjuanfuca[ , , 2] = DirStats::to_cir(data$zone_B)
sanjuanfuca[ , , 3] = DirStats::to_cir(data$zone_C)
sanjuanfuca[ , , 4] = DirStats::to_cir(data$zone_D)

# Save the object
save(list = "sanjuanfuca", file = "sanjuanfuca.rda")



