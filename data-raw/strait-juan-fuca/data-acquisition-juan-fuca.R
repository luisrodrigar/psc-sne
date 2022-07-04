
library(ncdf4)
library(lubridate)
library(circular)
library(dplyr)
library(here)
library(reshape2)
library(parallel)
library(abind)
library(tidyr)
library(ggmap)
library(viridis)

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
begin_lat <- 47.90
end_lat <- 48.8
begin_lon <- -124.20
end_lon <- -123.40

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


download_data <- TRUE # Change to TRUE download again in /data folder
if (download_data) {

  # Download data in a monthly loop to avoid surpassing the 500 Mb limit
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
        paste(loc, "_", toString(year), "_", toString(month), ".RData",
          sep = ""
        ),
        sep = "/"
      ))

      # Update month
      month <- month + 1

    }
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
  dir_results <- lapply(time_seq, function(i) {
    dir_speeds <- results[i:(i + day_ev - 1), c("d", "speed")]
    dir_speeds <- dir_speeds[complete.cases(dir_speeds), ]
    if (sum(dir_speeds$speed) == 0) {
      return(c(NA, NA))
    }
    weights <- dir_speeds$speed / sum(dir_speeds$speed)
    theta <- circular:::WeightedMeanCircularRad(w = weights, x = dir_speeds$d)
    speed <- mean(dir_speeds$speed)
    c(theta, speed)
  })
  return(do.call(rbind, dir_results))
}

# Remove rows by latitude and longitude with duplicate time entry
# Identified this occurs for the last and the first rows of consecutive months
remove_time_duplicates <- function(results_filter_latlon) {
  # Ids with the first duplicate row
  ids_to_remove <- results_filter_latlon %>%
    select(time) %>%
    duplicated()
  # Return the data.frame without this columns
  results_filter_latlon[!ids_to_remove, ]
}

# Convert to time instants of x hours where columns are location, time instant
# and theta
extract_long_fmt <- function(results, hours) {

  # Value for latitude and longitude
  latitudes <- levels(as.factor(results$lat))
  longitudes <- levels(as.factor(results$lon))

  # sequence from minimum to maximum time by each x hours
  instant <- seq(
    from = as.POSIXlt(min(results$time), tz = "UTC"),
    to = as.POSIXlt(max(results$time), tz = "UTC"),
    by = paste(hours, "hours", sep = " ")
  )

  # Grid of latitude and longitude
  grid <- expand.grid(latitudes, longitudes)

  # Calculate the weight mean circular theta each x hours
  location_results <- mclapply(
    mc.cores = detectCores() - 1, X = 1:nrow(grid),
    FUN = function(i) {
      # Extract latitude and longitude
      lat_x <- grid[i, 1]
      lon_x <- grid[i, 2]
      # Calculate the weight mean circular theta (x hours)
      theta_x <- results %>%
        filter(lat == lat_x & lon == lon_x) %>%
        remove_time_duplicates() %>%
        extract_theta(hours = hours)
      # Create data.frame
      data.frame(instant, lat_x, lon_x, theta_x)
    }
  )

  # Merge all the data.frame from above
  data <- do.call(rbind, location_results)
  # Arrange by ascending instant and rename column vars
  data %>%
    arrange(instant) %>%
    rename(
      time = instant,
      lat = lat_x,
      lon = lon_x,
      theta = X1,
      speed = X2
    )
}

# List individual RDatas
files <- list.files(
  path = here("data-raw", "strait-juan-fuca", "data"),
  pattern = "*.RData", full.names = TRUE, recursive = FALSE
)

# Obtain the results for every area and their corresponding daily directions
results <- extract_data(begin_lat, end_lat, begin_lon, end_lon)

# Data every 3 hours
hours <- 3
seq_time <- seq(
  from = as.POSIXlt(min(results$time), tz = "UTC"),
  to = as.POSIXlt(max(results$time), tz = "UTC"),
  by = paste(hours, "hours", sep = " ")
)

# Convert to theta by location (lat and lon) and time instant each 3 hours
juanfuca <- extract_long_fmt(results, hours = hours)

# Values of latitude and longitude
lon_values <- as.numeric(levels(juanfuca$lon))
lat_values <- as.numeric(levels(juanfuca$lat))

# 3-dimensional array to store the matrix of thetas for the latitude and longitude
lat_lon_theta_by_time <- array(
  dim = c(length(lon_values), length(lat_values), length(unique(juanfuca$time)))
)
# Select only the theta value with the location and time values
jdf_by_time <- juanfuca %>%
  select(time, lat, lon, theta)
# Create the 3-dimensional array split by time
jdf_by_time <- abind(split(jdf_by_time, jdf_by_time$time), along = 3)
for (k in 1:dim(jdf_by_time)[3]) {
  # Conver to wide format data frame
  # Where columns are longitude values and rows are latitude values
  lat_lon_theta_by_time[, , k] <- data.frame(jdf_by_time[, 2:4, k]) %>%
    reshape(timevar = "lon", idvar = "lat", direction = "wide") %>%
    select(-lat) %>%
    as.matrix()
}

# How many NA's are there?
filled.contour(lon_values, lat_values,
  apply(lat_lon_theta_by_time, 1:2, function(x) mean(is.na(x))),
  ylab = "lat", xlab = "lon", zlim = c(0, 1),
  main = "Proportion NAs", plot.axes = {
    box()
    axis(1)
    axis(2)
    points(expand.grid(lon_values, lat_values), pch = 16, cex = 0.5)
  }
)

# Calculate the percentage of NA's among the third dimension
percetange_na <- colMeans(apply(lat_lon_theta_by_time, is.na, MARGIN = c(1, 2)))
# Reflect the resulting matrix, reversing each row vector
percetange_na <- sapply(1:nrow(percetange_na), function(i) rev(percetange_na[i, ]))
# Adding rownames (latitude) and colnames (longitude)
rownames(percetange_na) <- rev(lat_values)
colnames(percetange_na) <- lon_values
# Convert the matrix lat x lon into long format
perc_na_long_fmt <- data.frame(percetange_na) %>%
  add_rownames(var = "lat") %>%
  data.frame %>%
  reshape(varying = (1:ncol(percetange_na) + 1), idvar = 'lat', timevar = "lon",
          direction = 'long', times = lon_values, v.names = "freq") %>%
  mutate(
    lat = round(as.numeric(lat), digits = 4),
    lon = round(lon, digits = 4)
  )

# Download the map, margin is increased
(map <- get_map(c(left = min(lon_values) - 0.05,
          bottom = min(lat_values) - 0.05,
          right = max(lon_values) + 0.05,
          top = max(lat_values) + 0.05)))

# Plot the map with the points and the countour surface colored
ggmap(map, extent = "panel") +
  ggtitle("Percentage of missing values") +
  geom_point(data = perc_na_long_fmt, aes(x = lon, y = lat, color = freq), size = 0.2) +
  geom_contour_filled(data = perc_na_long_fmt, aes(z = freq), alpha = 0.4, size = 0.1) +
  scale_color_viridis()


# Save the object
save(
  list = "juanfuca",
  file = paste(here("data-raw", "strait-juan-fuca"), "juanfuca.rda", sep = "/")
)

juanfuca_pkg <- paste(here("data-raw", "strait-juan-fuca"), "juanfuca.rda", sep = "/")
load(juanfuca_pkg)

# Transform to wide format the juan de fuca long data
juanfuca_wide <- juanfuca %>%
  mutate(
    location = paste(lat, lon, sep = ",")
  ) %>%
  select(-lat, -lon) %>%
  reshape(idvar = "time", timevar = "location", direction = "wide")

# Percentage of missing values
percetange_na_by_col <- colMeans(is.na(juanfuca_wide[, 2:ncol(juanfuca_wide)]))

# Set the maximum percentage of NA's, 40%
max_percentage_na <- 0.4
# Obtain only the columns that contains less than x% of NA's
jdf <- juanfuca_wide[, c(TRUE, percetange_na_by_col < max_percentage_na)]

# all non repeated locations from the wide data.frame
locations <- colnames(jdf)[seq(2, ncol(jdf), by = 2)]

# Obtaining the location names, use in reshape times parameter
time_names <-
  sapply(
    seq_along(locations),
    function(i) strsplit(locations[i], "theta.")[[1]][2]
  )

# long format with less than a percentage of NA's in location
jdf_long_fmt <- reshape(jdf,
  direction = "long", idvar = "time", timevar = "location",
  varying = 2:ncol(jdf), v.names = c("theta", "speed"), sep = ".",
  times = time_names
) %>%
  # separate location by comma separator
  separate(col = location, into = c("lat", "lon"), sep = ",") %>%
  rename(
    # reshape are putting incorrectly the order of each value for speed and theta
    # swapping values between speed and theta
    theta = speed,
    speed = theta
  ) %>%
  select(time, lat, lon, theta, speed)
# Remove index names
rownames(jdf_long_fmt) <- NULL

# Save dataset in long format
save(
  list = "jdf_long_fmt",
  file = paste(here("data-raw", "strait-juan-fuca"), "jdf.rda", sep = "/")
)
