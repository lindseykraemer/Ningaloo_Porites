#Climatological Temperature Range (sourced from NOAA's CRW) https://coralreefwatch.noaa.gov/product/thermal_history/climatology.php


#install required packages
#install.packages("ncdf4")

#################### Rough Draft on one site ########################


#load the libraries
library(ncdf4)

#set the path to downloaded netCDF file
file_path <- "noaa_crw_thermal_history_climatology_v3.5.nc"

#open the NetCDF file
nc_data <- nc_open(file_path)

#explore the data 

#list all variables in the file
print(nc_data)

#list variable names
names(nc_data$var)

#loop through variables and print their metadata
for (var_name in names(nc_data$var)) {
  cat(sprintf("\nVariable: %s\n", var_name))
  print(nc_data$var[[var_name]])
}

#extract data for 'clim_range' for specific location

#extract latitude and longitude arrays
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")

#extract the clim_range variable
clim_range <- ncvar_get(nc_data, "clim_range")

#target location (Bundegi)
target_lat <- -21.861278
target_lon <- 114.158750

#find the nearest indices for latitude and longitude
lat_index <- which.min(abs(lat - target_lat))
lon_index <- which.min(abs(lon - target_lon))

#get the closest latitude and longitude
closest_lat <- lat[lat_index]
closest_lon <- lon[lon_index]
cat(sprintf("Closest Latitude: %f, Closest Longitude: %f\n", closest_lat, closest_lon))

#extract the clim_range value
clim_range_value <- clim_range[lon_index, lat_index]
cat(sprintf("Climatological Range at (%f, %f): %.2f°C\n", closest_lat, closest_lon, clim_range_value))

#calculate SE in climatological temp range

#the clim_range is calculated as the difference between the maximum and minimum monthly mean SST climatology values (clim_monthly). To calculate SE:
#extract the monthly climatology values (clim_monthly)
#compute the range for all months at the specific location
#use the variability in these values to calculate the SE.

#extract the monthly SST climatology values
clim_monthly <- ncvar_get(nc_data, "clim_monthly")  #3D array [lon, lat, months]

#extract monthly values for the specific location
monthly_values <- clim_monthly[lon_index, lat_index, ]

#gut check: print the monthly SST values
cat("Monthly SST values (°C):\n")
print(monthly_values)

#calculate the climatological range
clim_range_calc <- max(monthly_values, na.rm = TRUE) - min(monthly_values, na.rm = TRUE)

#calculate the SD of monthly SST
std_dev <- sd(monthly_values, na.rm = TRUE)

#calculate the SE
se <- std_dev / sqrt(length(monthly_values))

#print results
cat(sprintf("Climatological Range (calculated): %.2f°C\n", clim_range_calc))
cat(sprintf("Standard Error (SE): %.2f°C\n", se))

#compare calculated clim_range with the value provided in the NetCDF file to verify consistency.

#compare with clim_range from the file
cat(sprintf("Climatological Range (NetCDF): %.2f°C\n", clim_range_value))

#good to go

############# Now for all sites #############


#define the locations with their *decimal degree* coordinates
locations <- data.frame(
  Name = c("BUN", "VLF", "OSP", "TAN", "MUI", "COR", "LAK/TUR"),
  Latitude = c(-21.861278, -21.791944, -22.242778, -21.899167, -21.677889, -23.135694, -22.06750),
  Longitude = c(114.15875, 114.17500, 113.828889, 113.968333, 114.341361, 113.76400, 113.897222)
)

#extract and analyse data for each location

#load libraries
library(ncdf4)

#open the NetCDF file
file_path <- "noaa_crw_thermal_history_climatology_v3.5.nc"
nc_data <- nc_open(file_path)

#extract latitude, longitude, and clim_monthly variable
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")
clim_range <- ncvar_get(nc_data, "clim_range")
clim_monthly <- ncvar_get(nc_data, "clim_monthly")

#initialize results df
results <- data.frame(
  Name = character(),
  Closest_Lat = numeric(),
  Closest_Lon = numeric(),
  Clim_Range = numeric(),
  SE = numeric(),
  stringsAsFactors = FALSE
)

#loop through each location
for (i in 1:nrow(locations)) {
  # Get target coordinates
  target_lat <- locations$Latitude[i]
  target_lon <- locations$Longitude[i]
  
  #find the nearest indices
  lat_index <- which.min(abs(lat - target_lat))
  lon_index <- which.min(abs(lon - target_lon))
  
  #get the closest grid coordinates
  closest_lat <- lat[lat_index]
  closest_lon <- lon[lon_index]
  
  #extract clim_range value
  clim_range_value <- clim_range[lon_index, lat_index]
  
  #extract monthly climatology values
  monthly_values <- clim_monthly[lon_index, lat_index, ]
  
  #calculate SE
  std_dev <- sd(monthly_values, na.rm = TRUE)
  se <- std_dev / sqrt(length(monthly_values))
  
  #results
  results <- rbind(results, data.frame(
    Name = locations$Name[i],
    Closest_Lat = closest_lat,
    Closest_Lon = closest_lon,
    Clim_Range = clim_range_value,
    SE = se
  ))
}

#print results
print(results)

#save results to a CSV file
write.csv(results, "clim_range_results.csv", row.names = FALSE)

#close the NetCDF file
nc_close(nc_data)
