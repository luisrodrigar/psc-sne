
library(ncdf4)
library(lubridate)
library(circular)
library(dplyr)
library(here)
library(reshape2)

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
begin_lat <- 48.175
end_lat <- 48.392
begin_lon <- -124.18
end_lon <- -123.48

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
    units = "Mb"
  ))

  # Download (u, v), being:
  # u (m/s) = surface_eastward_sea_water_velocity
  # v (m/s) = surface_northward_sea_water_velocity
  u <<- ncvar_get(data, "u",
    start = c(begin_lon_ind, begin_lat_ind, begin_tim_ind),
    count = c(l_lon, l_lat, l_tim)
  )
  v <<- ncvar_get(data, "v",
    start = c(begin_lon_ind, begin_lat_ind, begin_tim_ind),
    count = c(l_lon, l_lat, l_tim)
  )
}

# Download data in a monthly loop to avoid surpassing 500Mb limit
for (year in 2020:2022) {
  month <- 1
  while ((year < 2022 && month <= 12) || (year == 2022 && month <= 6)) {

    # Show progress
    cat("\n", year, "-", month, "\n")

    begin_tim <- date(paste(toString(year), toString(month),
      "01 00:00:00 UTC",
      sep = "-"
    ))

    # End time taking into account if a new year starts
    if (month != 12) {
      end_tim <- date(paste(toString(year), toString(month + 1),
        "01 00:00:00 UTC",
        sep = "-"
      ))
    } else {
      end_tim <- date(paste(toString(year + 1), toString(1),
        "01 00:00:00 UTC",
        sep = "-"
      ))
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
    save(final_dataframe, file = paste(
      here("data-raw", "strait-juan-fuca", "data"),
      paste(loc, "_", toString(year), "_", toString(month), ".RData", sep = ""),
      sep = "/"
    ))

    # Update month
    month <- month + 1
  }
}

## Obtaining daily data for the specific zone

# Function that reads over all the files in the directory containing the raw
# data and return the records inside a given area delimited by longitude
# and latitude
extract_data <- function(begin_lat, end_lat, begin_lon, end_lon) {

  # Retrieve monthly data
  monthly_data <- lapply(files, function(x) {

    # Load raw data file
    load(x)
    results <- filter(
      final_dataframe, lat > begin_lat, lat < end_lat,
      lon > begin_lon, lon < end_lon
    )

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

  # Dimensions
  n <- length(levels(as.factor(results$lat)))
  m <- length(levels(as.factor(results$lon)))
  t <- hours
  day_ev <- n * m * t
  total_events <- length(results$u)


  # Apply the speed-weighted average over the non-NA observations inside the
  # t hour period
  time_seq <- seq(1, total_events, day_ev)
  theta <- sapply(time_seq, function(i) {
    dir_speeds <- results[i:(i + day_ev - 1), c("d", "speed")]
    dir_speeds <- dir_speeds[complete.cases(dir_speeds), ]
    if(sum(dir_speeds$speed) == 0) {
      return(NA)
    }
    weights <- dir_speeds$speed / sum(dir_speeds$speed)
    circular:::WeightedMeanCircularRad(w = weights, x = dir_speeds$d)
  })
  return(theta)
}

# Remove rows by latitude and longitude with duplicate time entry
# Identified this occurs for the last and the first rows of consecutive months
remove_time_duplicates <- function(results_filter_latlon) {
  # Ids with the first duplicate row
  ids_to_remove <- results_filter_latlon %>%
    select(time) %>%
    duplicated
  # Return the data.frame without this columns
  results_filter_latlon[!ids_to_remove, ]
}

# Convert to time instants of x hours where columns are location, time instan and theta
theta_by_inst_location <- function(results, hours) {
  # Value for latitude and longitude
  latitudes <- levels(as.factor(results$lat))
  longitudes <- levels(as.factor(results$lon))

  # sequence from minimum to maximum time by each x hours
  instant <- seq(from = as.POSIXlt(min(results$time), tz="UTC"),
                       to = as.POSIXlt(max(results$time), tz="UTC"),
                       by = paste(hours, 'hours', sep = " "))
  # Grid of latitude and longitude
  grid <- expand.grid(latitudes, longitudes)
  # Calculate the weight mean circular theta each x hours
  location_results <- mclapply(mc.cores = detectCores() - 1, X = 1:nrow(grid),
                               FUN = function(i) {
                                 # Extract latitude and longitude
                                 lat_x <- grid[i, 1]
                                 lon_x <- grid[i, 2]
                                 # Create location based on the two previous values
                                 location_x <- paste(lat_x, lon_x, sep = ",")
                                 # repeat y times the location, where y is the size of the time interval
                                 location_x <- rep(location_x, size)
                                 # Calculate the weight mean circular theta (x hours)
                                 theta_x <- results %>%
                                   filter(lat == lat_x & lon == lon_x) %>%
                                   remove_time_duplicates %>%
                                   extract_theta(., hours = hours)
                                 # Create data.frame
                                 data.frame(instant, location_x, theta_x)
                              })
  # Merge all the data.frame from above
  data <- do.call(rbind, location_results)
  return(data)
}

# List individual RDatas
files <- list.files(path = here("data-raw", "strait-juan-fuca", "data"),
                    pattern = "*.RData", full.names = TRUE, recursive = FALSE)

# Obtain the results for every area and their corresponding daily directions
results <- extract_data(begin_lat, end_lat, begin_lon, end_lon)
# Data every 3 hours
hours <- 3
# Theta by location and time instant each 3 hours
theta_inst_location <- theta_by_inst_location(results = results, hours = hours)
# Convert to matrix where each row is a time instant.
# Each column is a different location (latitude and longitude)
juan_fuca <- reshape2::dcast(data = theta_inst_location,
                   formula = instant ~ location_x,
                   value.var = 'theta_x', fun.aggregate = list)
rownames(juan_fuca) <- seq(from = as.POSIXlt(min(results$time), tz="UTC"),
                           to = as.POSIXlt(max(results$time), tz="UTC"),
                           by = paste(hours, 'hours', sep = " "))


# Save the object
save(list = "juan_fuca",
     file = paste(here('data-raw', 'strait-juan-fuca'), "juanfuca.rda", sep = "/"))
