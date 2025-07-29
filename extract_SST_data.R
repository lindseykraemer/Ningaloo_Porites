# Code for determining historic mean SST at each site
# Also determining which site has historically been the warmest or coolest on average

#load libraries
library(ncdf4)

#path to the climatology dataset
file_path <- "noaa_crw_thermal_history_climatology_v3.5.nc"

#open the NetCDF file
nc_data <- nc_open(file_path)

#explore data structure

#list all variables in the file
print(nc_data)

#list variable names
names(nc_data$var)

#loop through variables and print their metadata
for (var_name in names(nc_data$var)) {
  cat(sprintf("\nVariable: %s\n", var_name))
  print(nc_data$var[[var_name]])
}


#check available years in database
years <- ncvar_get(nc_data, "years")  #extract years

#print the extracted time period
print(years)  #list of years analysed

#extract global attributes related to time coverage
time_start <- ncatt_get(nc_data, 0, "time_coverage_start")$value
time_end <- ncatt_get(nc_data, 0, "time_coverage_end")$value

#print the time period analyzed
cat("Time period analyzed:", time_start, "to", time_end, "\n")

#extract latitude and longitude grids
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")

#extract monthly climatology SST variable
clim_monthly <- ncvar_get(nc_data, "clim_monthly")  # [lon, lat, months]

#define the locations with their decimal degree coordinates
locations <- data.frame(
  Name = c("BUN", "VLF", "OSP", "TAN", "MUI", "COR", "LAK/TUR"),
  Latitude = c(-21.861278, -21.791944, -22.242778, -21.899167, -21.677889, -23.135694, -22.06750),
  Longitude = c(114.15875, 114.17500, 113.828889, 113.968333, 114.341361, 113.76400, 113.897222)
)

#extract data

#initialize results data frame
results <- data.frame(
  Name = character(),
  Closest_Lat = numeric(),
  Closest_Lon = numeric(),
  Annual_Mean_SST = numeric(),
  SE_Annual_Mean_SST = numeric(),
  Warmest_Month = numeric(),
  Warmest_Month_SST = numeric(),
  SE_Warmest_Month_SST = numeric(),
  Coolest_Month = numeric(),
  Coolest_Month_SST = numeric(),
  SE_Coolest_Month_SST = numeric(),
  stringsAsFactors = FALSE
)

#loop through each location
for (i in 1:nrow(locations)) {
  #get target location
  target_lat <- locations$Latitude[i]
  target_lon <- locations$Longitude[i]
  
  #find the nearest grid point
  lat_index <- which.min(abs(lat - target_lat))
  lon_index <- which.min(abs(lon - target_lon))
  
  #extract SST values for the location
  sst_values <- clim_monthly[lon_index, lat_index, ]  #monthly SST for the site
  
  #calculate metrics
  annual_mean_sst <- mean(sst_values, na.rm = TRUE)  #annual mean
  annual_sst_sd <- sd(sst_values, na.rm = TRUE)      #SD
  se_annual_mean_sst <- annual_sst_sd / sqrt(12)     #SE
  
  warmest_month <- which.max(sst_values)             #warmest month (1–12)
  coolest_month <- which.min(sst_values)             #coolest month (1–12)
  
  warmest_month_sst <- sst_values[warmest_month]
  coolest_month_sst <- sst_values[coolest_month]
  
  
  se_warmest_month_sst <- annual_sst_sd / sqrt(1)  #1 observation per warmest month
  se_coolest_month_sst <- annual_sst_sd / sqrt(1)  #1 observation per coolest month
  
  
  #create a df with results
  results <- rbind(results, data.frame(
    Name = locations$Name[i],
    Closest_Lat = lat[lat_index],
    Closest_Lon = lon[lon_index],
    Annual_Mean_SST = annual_mean_sst,
    SE_Annual_Mean_SST = se_annual_mean_sst,
    Warmest_Month = warmest_month,
    Warmest_Month_SST = warmest_month_sst,
    SE_Warmest_Month_SST = se_warmest_month_sst,
    Coolest_Month = coolest_month,
    Coolest_Month_SST = coolest_month_sst,
    SE_Coolest_Month_SST = se_coolest_month_sst
  ))
}


print(results)


#save results
write.csv(results, "sst_analysis_results.csv", row.names = FALSE)

#slose the file
nc_close(nc_data)
