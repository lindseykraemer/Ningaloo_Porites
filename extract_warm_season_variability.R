#Code for warm season variability (sourced from NOAA CRW) https://coralreefwatch.noaa.gov/product/thermal_history/sst_variability.php

#noaa_crw_thermal_history_sst_variability_v3.5

################### Now for all sites ###################

#including 'mean_withinyear_sd_warmseason' provides the standard deviation of warm season variability for each year
#the number of years (n) is derived from the file metadata (39 years for 1985–2023)


#define locations in a data frame
locations <- data.frame(
  Name = c("BUN", "VLF", "OSP", "TAN", "MUI", "COR", "LAK/TUR"),
  Latitude = c(-21.861278, -21.791944, -22.242778, -21.899167, -21.677889, -23.135694, -22.06750),
  Longitude = c(114.15875, 114.17500, 113.828889, 113.968333, 114.341361, 113.76400, 113.897222)
)

#load libraries
library(ncdf4)

#open the NetCDF file
file_path <- "noaa_crw_thermal_history_sst_variability_v3.5.nc"
nc_data <- nc_open(file_path)

#extract latitude, longitude, and warm season variability data
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")
warm_season_var <- ncvar_get(nc_data, "stdv_warmseason")
mean_withinyear_sd_warmseason <- ncvar_get(nc_data, "mean_withinyear_sd_warmseason")  # Optional for SE

#initialize results data frame
results <- data.frame(
  Name = character(),
  Closest_Lat = numeric(),
  Closest_Lon = numeric(),
  Warm_Season_Var = numeric(),
  SE_Warm_Season_Var = numeric(),
  stringsAsFactors = FALSE
)

#loop through each location
for (i in 1:nrow(locations)) {
  #get target coordinates
  target_lat <- locations$Latitude[i]
  target_lon <- locations$Longitude[i]
  
  #find nearest grid point
  lat_index <- which.min(abs(lat - target_lat))
  lon_index <- which.min(abs(lon - target_lon))
  
  #extract warm season variability and calculate SE
  warm_season_value <- warm_season_var[lon_index, lat_index]
  
  #calculate SE based on mean_withinyear_sd_warmseason 
  se_warm_season <- mean_withinyear_sd_warmseason[lon_index, lat_index] / sqrt(39)  #divide by sqrt of number of years
  
  #append results
  results <- rbind(results, data.frame(
    Name = locations$Name[i],
    Closest_Lat = lat[lat_index],
    Closest_Lon = lon[lon_index],
    Warm_Season_Var = warm_season_value,
    SE_Warm_Season_Var = se_warm_season
  ))
}

#print results
print(results)

#save results to CSV
write.csv(results, "warm_season_variability.csv", row.names = FALSE)

#close the NetCDF file
nc_close(nc_data)

