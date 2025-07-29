# warm season variability site map (bounding box)


#load libraries and open NetCDF file
library(ncdf4)
library(dplyr)

#open the NetCDF file
file_path <- "noaa_crw_thermal_history_sst_variability_v3.5.nc"
nc_data <- nc_open(file_path)

#extract latitude, longitude, relevant variables, and the mask
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")
warm_season_var <- ncvar_get(nc_data, "stdv_warmseason")
mean_withinyear_sd_warmseason <- ncvar_get(nc_data, "mean_withinyear_sd_warmseason")  #for SE
mask <- ncvar_get(nc_data, "mask")  #use this to exclude non-analyzed pixels

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
  Warm_Season_Var = numeric(),
  SE = numeric(),
  stringsAsFactors = FALSE
)

#loop through all grid cells within the bounding box
for (i in lon_indices) {
  for (j in lat_indices) {
    #check if the grid point is valid based on the mask
    if (mask[i, j] == 0) next  #skip non-analyzed pixels
    
    #get the closest latitude and longitude
    closest_lon <- lon[i]
    closest_lat <- lat[j]
    
    #extract warm season variability value
    warm_season_value <- warm_season_var[i, j]
    
    #skip NA values or values below the valid minimum
    if (is.na(warm_season_value) || warm_season_value < 0.169626474380493) next
    
    #calculate SE based on mean_withinyear_sd_warmseason
    se <- mean_withinyear_sd_warmseason[i, j] / sqrt(39)  # Assuming 39 years of data
    
    #append to results
    results <- rbind(results, data.frame(
      Lon = closest_lon,
      Lat = closest_lat,
      Warm_Season_Var = warm_season_value,
      SE = se
    ))
  }
}

#print results
print(head(results))

#save results to a CSV file
write.csv(results, "warm_season_results_bounding_box.csv", row.names = FALSE)

#close the NetCDF file
nc_close(nc_data)

#made sure this aligned with same data for specific locations from warm_seas_var data folder
