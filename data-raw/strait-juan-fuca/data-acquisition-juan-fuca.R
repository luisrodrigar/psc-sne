
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
begin_lat <- 47.85
end_lat <- 48.64
begin_lon <- -124.254
end_lon <- -123.36

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
  # objects have dimension c(l_lon, l_lat, l_tim)
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
    if(year == 2020) {

      # Starting in June for the year 2020
      month <- 6

    }

    while ((year <= 2021 && month <= 12) || (year == 2022 && month <= 6)) {

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
      lat_length <- dim(u)[2]
      lon_length <- dim(u)[1]
      time_length <- dim(u)[3]

      # Find the indexes associated to dimensions
      lat_aux <- begin_lat_ind:(begin_lat_ind + lat_length - 1)
      lon_aux <- begin_lon_ind:(begin_lon_ind + lon_length - 1)
      time_aux <- begin_tim_ind:(begin_tim_ind + time_length - 1)

      # Join all the cases
      # Since u and v are objects of dimension c(l_lon, l_lat, l_tim)
      # The combinations must be done between each latitude and its possible
      # combinatin with the longitude values. Later on, those three-dimensional
      # arrays will be passed as vector (c(array) converts into a vector first
      # by each column, later each row and finally the third dimension)
      join <- merge(x = lon[lon_aux], y = lat[lat_aux])
      final_dataframe <- merge(x = join, y = tim[time_aux], by = NULL)
      colnames(final_dataframe) <- c("lon", "lat", "time")

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

# Remove rows by latitude and longitude with duplicate time entry
# Identified this occurs for the last and the first rows of consecutive months
remove_next_month_instants <- function(results) {

  # Get the set of time instants
  dates <- unique(as.POSIXct(results$time, tz = "UTC"))
  # Convert them into their associated month
  months <- month(dates)
  # Calculate the frequency table
  freq_table <- table(months)

  # Only if there is more than one month in the data
  if (length(freq_table) > 1) {

    # Min and max dates in a unit dataset which is a month interval
    least_frequent_index <- which(freq_table == min(freq_table))
    least_frequent_month <- as.numeric(names(freq_table[least_frequent_index]))

    # Get dates to remove from original dataset
    date_to_remove <- dates[which(months == least_frequent_month)]
    # Exclude that date
    results <- results %>%
      filter(!(as.POSIXct(time, tz = "UTC") %in% as.POSIXct(date_to_remove, tz = "UTC")))

  }
  return(results)
}


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

    # Remove next month instants
    results <- remove_next_month_instants(results)

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
extract_theta <- function(results, hours, lat, lon, instant) {

  # Dimensions
  t <- hours
  total_events <- length(results$u)

  # Apply the speed-weighted average over the non-NA observations inside the
  # t hour period
  time_seq <- seq(1, total_events, t)
  data_t_time <- lapply(
    time_seq,
    function(i) {
      # Obtaining the rows within the t interval of time
      dir_speeds <- results[i:(i + t - 1), c("d", "speed")]
      # Removing rows that contain any NA's
      dir_speeds <- dir_speeds[complete.cases(dir_speeds), ]
      # The index of the time instant
      instant_index <- ceiling(i / hours)
      # return NA for empty dir_speed or sum of speeds equal to 0
      if (is.null(dir_speeds$speed) || sum(dir_speeds$speed) == 0) {

       # Returning the values NA for the direction angle and speed
       return(c(as.numeric(instant[instant_index]), lat, lon, NA, NA))

      } else {

       # Return the weight mean circular direction
       # over the t interval and the speed mean
       return(c(as.numeric(instant[instant_index]), lat, lon,
                circular:::WeightedMeanCircularRad(
                  w = dir_speeds$speed / sum(dir_speeds$speed),
                  x = dir_speeds$d),
                mean(dir_speeds$speed)))

      }
    })
  return(do.call(rbind, data_t_time))
}


# Convert to time instants of x hours where columns are location, time instant
# and theta
extract_long_fmt <- function(results, hours, num_cores = detectCores() - 1) {

  # sequence from minimum to maximum time by each x hours
  instant_seq <- seq(
    from = as.POSIXlt(min(results$time), tz = "UTC"),
    to = as.POSIXlt(max(results$time), tz = "UTC"),
    by = paste(hours, "hours", sep = " ")
  )

  # Grid of latitude and longitude
  grid <- expand.grid(unique(results$lat), unique(results$lon))

  location_results <- mclapply(
    mc.cores = num_cores,
    X = 1:nrow(grid),
    FUN = function(i) {
      results %>%
        filter(lat == grid[i, 1] & lon == grid[i, 2]) %>%
        extract_theta(hours = hours,
                      lat = grid[i, 1],
                      lon = grid[i, 2],
                      instant = instant_seq)
      })

  # Merge all the data.frame from above
  data <- do.call(rbind, location_results)
  # Arrange by ascending instant and rename column vars
  data.frame(data) %>%
    mutate(
      time = as.POSIXct(.[[1]], origin = '1970-01-01', tz = 'UTC'),
      lat = .[[2]],
      lon = .[[3]],
      theta = .[[4]],
      speed = .[[5]]
    ) %>%
    dplyr::select(time, lat, lon, theta, speed) %>%
    arrange(time)
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

# Convert to theta by location (lat and lon) and time instant each 3 hours
juanfuca <- extract_long_fmt(results, hours = hours, num_cores = 7)

# Save the object
save(
  list = "juanfuca",
  file = paste(here("data-raw", "strait-juan-fuca"), "juanfuca.rda", sep = "/")
)

# Remove source data
rm(results)
# Garbage collector: clean unused memory space
gc()

load(paste(here("data-raw", "strait-juan-fuca"), "juanfuca.rda", sep = "/"))

theta_lat_lon_create <- function(lon_values, lat_values, data) {

  # Select only the theta value with the location and time values
  jdf_by_time <- data %>%
    dplyr::select(time, lat, lon, theta) %>%
    mutate(time = as.character(time, format="%Y-%m-%d %H:%M:%S"),
           lat = as.character(lat),
           lon = as.character(lon),
           theta = as.character(theta))
  # Create the 3-dimensional array split by time
  jdf_by_time <- abind(split(jdf_by_time, jdf_by_time$time), along = 3)

  # Time instant size
  r <- dim(jdf_by_time)[3]

  # 3-dimensional array to store the matrix of thetas for the latitude and longitude
  lat_lon_theta_by_time <- array(
    dim = c(length(lat_values), length(lon_values), r)
  )

  for (k in seq_len(r)) {
    # Convert to wide format data frame
    # Where columns are longitude values and rows are latitude values
    lat_lon_theta_by_time[, , k] = (data.frame(jdf_by_time[ , 2:4, k]) %>%
                                      mutate(lat = as.numeric(lat),
                                             lon = as.numeric(lon),
                                             theta = as.numeric(theta)) %>%
                                      reshape(timevar = "lon", idvar = "lat", direction = "wide") %>%
                                      dplyr::select(-lat) %>%
                                      as.matrix())
  }
  return(lat_lon_theta_by_time)
}

get_percentage_na <- function(data) {

  return(data %>%
           group_by(lat, lon) %>%
           summarise(freq = mean(is.na(theta))) %>%
           data.frame)
}

# Values of latitude and longitude
lon_values <- unique(juanfuca$lon)
lat_values <- unique(juanfuca$lat)

# Create objet to use in filled.contour
lat_lon_theta_by_time <- theta_lat_lon_create(lon_values = lon_values,
                                              lat_values = lat_values,
                                              data = juanfuca)

# How many NA's are there?
filled.contour(lon_values, lat_values,
               t(apply(lat_lon_theta_by_time, 1:2, function(x) mean(is.na(x)))),
               ylab = "lat", xlab = "lon", zlim = c(0, 1),
               main = "Proportion NAs", plot.axes = {
                 box()
                 axis(1)
                 axis(2)
                 points(expand.grid(lon_values, lat_values), pch = 16, cex = 0.5)
                 })

# Show percentage of missing values on the map (Strait Juan de Fuca)
perc_na_long_fmt <- get_percentage_na(juanfuca)

# Download the map, margin is increased
(map <- get_map(c(left = min(lon_values) - 0.05,
                  bottom = min(lat_values) - 0.05,
                  right = max(lon_values) + 0.05,
                  top = max(lat_values) + 0.05)))

# Plot the map with the points and the countour surface colored
ggmap(map, extent = "panel") +
  ggtitle("Percentage of missing values") +
  geom_contour_filled(data = perc_na_long_fmt, aes(z = freq), alpha = 0.6, size = 0.1)

# Merge the previous dataset that contains the frequency of missing values
# with the base dataset that contains all the information for this interval of time
jdf_freq <- merge(x = juanfuca, y = perc_na_long_fmt, by = c("lat", "lon"))

# Set the maximum percentage of NA's, 10%
max_percentage_na <- 0.10

jdf_long_fmt <- jdf_freq %>%
  filter(freq <= max_percentage_na) %>%
  arrange(time)

# Plot in the map the locations selected
ggmap(map, extent = "panel") +
  ggtitle("Selected data locations for the analysis") +
  geom_contour_filled(data = jdf_long_fmt, aes(z = freq), alpha = 0.6, size = 0.1)

# Remove freq column
jdf_long_fmt <- jdf_long_fmt %>%
  select(-freq)

# Save dataset in long format
jdf <- jdf_long_fmt
save(
  list = "jdf",
  file = paste(here("data-raw", "strait-juan-fuca"), "jdf.rda", sep = "/"),
  compress = "bzip2"
)
