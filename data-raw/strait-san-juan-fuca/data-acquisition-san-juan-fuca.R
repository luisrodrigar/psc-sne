
library(ncdf4)
library(lubridate)
library(circular)
library(dplyr)

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
for (year in 2020:2022) {

  month <- 1
  if (year == 2020) {
    month <- 5
  }
  while ( ( year < 2022 && month <= 12 ) || (year == 2022 && month <= 5) ) {

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

# List individual RDatas
files <- list.files(pattern = "*.RData", full.names = TRUE, recursive = FALSE)

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
extract_theta <- function(results, hours) {
  
  n <- length(levels(as.factor(results$lat)))
  m <- length(levels(as.factor(results$lon)))
  t <- hours
  day_ev <- n * m * t
  total_events <- length(results$u)

  ttt <- results %>% 
    mutate(
      month = lubridate::month(time), 
      year = lubridate::year(time),
      date = paste(year, month, sep="-"),
      location = paste(lat, lon, sep=",")
    ) %>% 
    group_by(lat, lon, year, month, .drop=T) %>% 
    na.omit %>% 
    mutate(
      weights = speed / sum(speed), 
      theta=circular:::WeightedMeanCircularRad(w = weights, x = d)
    ) %>% 
    ungroup %>% 
    select(location, date, theta) %>% 
    unique %>% 
    mutate(
      coord = DirStats::to_cir(theta)
    ) %>% 
    select(-theta) %>% 
    data.frame
  ttt <- ttt %>% filter(date != "2022-6")
  sanjuanfuca <- abind(split(ttt, ttt$date), along=3)
  rownames(sanjuanfuca) <- sanjuanfuca[,1,1]
  sanjuanfuca <- sanjuanfuca[,3:4,]

  # Apply the speed-weighted average over the non-NA observations inside the
  # t hour period
  time_seq <- seq(1, total_events, day_ev)
  results <- results[order(as.Date(results$time)),]
  day_seq <- as.character(seq(min(results$time), max(results$time), 
                              by=paste(hours, "hours")))
  theta <- sapply(time_seq, function(i) {
    
    dir_speeds <- results[i:(i + day_ev - 1), c("d", "speed")]
    dir_speeds <- dir_speeds[complete.cases(dir_speeds), ]
    weights <- dir_speeds$speed / sum(dir_speeds$speed)
    list(date=day_seq[as.integer(i / day_ev) + 1],
         theta=circular:::WeightedMeanCircularRad(w = weights, x = dir_speeds$d))
    
  }, simplify = 'array')
  return(t(theta))
  
}

results <- extract_data(begin_lat, end_lat, begin_lon, end_lon)
results <- results[order(as.Date(results$time)),]
ttt <- results %>% 
  mutate(
    month = lubridate::month(time), 
    year = lubridate::year(time),
    date = paste(year, month, sep="-"),
    location = paste(lat, lon, sep=",")
  ) %>% 
  group_by(lat, lon, year, month, .drop=T) %>% 
  na.omit %>% 
  mutate(
    weights = speed / sum(speed), 
    theta=circular:::WeightedMeanCircularRad(w = weights, x = d)
  ) %>% 
  ungroup %>% 
  select(location, date, theta) %>% 
  unique %>% 
  mutate(
    coord = DirStats::to_cir(theta)
  ) %>% 
  select(-theta) %>% 
  data.frame
ttt <- ttt %>% filter(date != "2022-6")
sanjuanfuca <- abind(split(ttt, ttt$date), along=3)
rownames(sanjuanfuca) <- sanjuanfuca[,1,1]
sanjuanfuca <- sanjuanfuca[,3:4,]
for (k in dim(sanjuanfuca)[3]) {
  sanjuanfuca[, 1:2, k] <- as.numeric(sanjuanfuca[, 1:2, k])
}

# Save the object
save(list = "sanjuanfuca", file = "sanjuanfuca.rda", compress = "bzip2")



