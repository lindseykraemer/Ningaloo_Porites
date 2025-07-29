#Code for determining NUMBER of Severe Bleaching-level Heat Stress Events (DHW≥8) and Time between Severe Bleaching-level Heat Stress Events (DHW≥8)
#(extracted from https://coralreefwatch.noaa.gov/product/thermal_history/stress_frequency.php)


#load libraries
library(ncdf4)
library(tidyverse)

#path to the climatology dataset
file_path <- "noaa_crw_thermal_history_stress_freq_v3.5.nc"

#open the NetCDF file
nc_data <- nc_open(file_path)

#list all variables in the file
print(nc_data)

#list variable names
names(nc_data$var)

#loop through variables and print their metadata
for (var_name in names(nc_data$var)) {
  cat(sprintf("\nVariable: %s\n", var_name))
  print(nc_data$var[[var_name]])
}

#"n_ge8" is the variable of interest for DHW > 8
#"rp_ge8" is the variable of interest for recovery time


#define variables

#extract variables
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")
rp_ge8 <- ncvar_get(nc_data, "rp_ge8")  #mean time between events
n_ge8 <- ncvar_get(nc_data, "n_ge8")    #number of events


#define the locations with their decimal degree coordinates
locations <- data.frame(
  Name = c("BUN", "VLF", "OSP", "TAN", "MUI", "COR", "LAK/TUR"),
  Latitude = c(-21.861278, -21.791944, -22.242778, -21.899167, -21.677889, -23.135694, -22.06750),
  Longitude = c(114.15875, 114.17500, 113.828889, 113.968333, 114.341361, 113.76400, 113.897222)
)

#find nearest grid points
find_nearest_index <- function(value, grid) {
  which.min(abs(grid - value))
}

locations <- locations %>%
  mutate(
    lon_index = map_dbl(Longitude, find_nearest_index, lon),
    lat_index = map_dbl(Latitude, find_nearest_index, lat)
  )


#extract metrics 
results <- locations %>%
  mutate(
    Time_Between_Severe_Events = map2_dbl(lon_index, lat_index, ~ rp_ge8[.x, .y]),
    Number_of_Events = map2_dbl(lon_index, lat_index, ~ n_ge8[.x, .y]),
    SD_Between_Events = Time_Between_Severe_Events * 0.15,  #assuming SD as 15% of mean
    SE_Between_Events = ifelse(Number_of_Events > 0, SD_Between_Events / sqrt(Number_of_Events), NA)
  )

#save results
write.csv(results, "DHW_events.csv", row.names = FALSE)
