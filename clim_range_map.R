# Gathering climatological temp range data for site map (bounding box)

# For the general map area, there is climatological temp range data missing for certain general areas. Incorporating interpolated data to account for this. 


#load libraries and open NetCDF file
library(ncdf4)
library(dplyr)

file_path <- "noaa_crw_thermal_history_climatology_v3.5.nc"
nc_data <- nc_open(file_path)

#extract latitude, longitude, climatological range, monthly climatology, and mask
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")
clim_range <- ncvar_get(nc_data, "clim_range")
clim_monthly <- ncvar_get(nc_data, "clim_monthly")  # Dimensions: [lon, lat, months]
mask <- ncvar_get(nc_data, "mask")  # Extract the mask variable

#define the latitudinal and longitudinal boundaries of the area
min_lon <- 113.0
max_lon <- 115.0
min_lat <- -24.0
max_lat <- -21.0

#find indices within the boundaries
lon_indices <- which(lon >= min_lon & lon <= max_lon)
lat_indices <- which(lat >= min_lat & lat <= max_lat)

#subset the longitude and latitude arrays
lon_subset <- lon[lon_indices]
lat_subset <- lat[lat_indices]

#initialize results dataframe
results <- data.frame(
  Lon = numeric(),
  Lat = numeric(),
  Clim_Range = numeric(),
  SE = numeric(),
  stringsAsFactors = FALSE
)

#loop through all grid cells within the bounding box
for (i in lon_indices) {
  for (j in lat_indices) {
    #skip non-analyzed pixels based on the mask
    if (mask[i, j] == 0) next
    
    #get the closest latitude and longitude
    closest_lon <- lon[i]
    closest_lat <- lat[j]
    
    #extract clim_range value
    clim_range_value <- clim_range[i, j]
    
    #extract monthly climatology values
    monthly_values <- clim_monthly[i, j, ]
    
    #skip NA values
    if (all(is.na(monthly_values))) next
    
    #calculate SE
    std_dev <- sd(monthly_values, na.rm = TRUE)
    se <- std_dev / sqrt(length(monthly_values))
    
    #append to results
    results <- rbind(results, data.frame(
      Lon = closest_lon,
      Lat = closest_lat,
      Clim_Range = clim_range_value,
      SE = se
    ))
  }
}

#print results
print(results)

#save the results to a CSV file
write.csv(results, "clim_range_results_bounding_box.csv", row.names = FALSE)

#close the NetCDF file
nc_close(nc_data)

# Made sure this aligned with same data for specific locations from climatological range data folder
