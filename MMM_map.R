# Code for creating a map based on MMM (Warmest_Month_SST)

#load libraries
library(ncdf4)
library(dplyr)

#open the NetCDF file
file_path <- "noaa_crw_thermal_history_climatology_v3.5.nc"
nc_data <- nc_open(file_path)

#extract latitude, longitude, mask, monthly SST, and warmest month
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")
clim_monthly <- ncvar_get(nc_data, "clim_monthly")  #[lon, lat, months]
mmm_month <- ncvar_get(nc_data, "mmm_month")        #[lon, lat]
mask <- ncvar_get(nc_data, "mask")                 #[lon, lat]

#define the bounding box
min_lon <- 113.0
max_lon <- 115.0
min_lat <- -24.0
max_lat <- -21.0

#find indices within the bounding box
lon_indices <- which(lon >= min_lon & lon <= max_lon)
lat_indices <- which(lat >= min_lat & lat <= max_lat)

#subset longitude and latitude
lon_subset <- lon[lon_indices]
lat_subset <- lat[lat_indices]

#initialize results data frame
results <- data.frame(
  Lon = numeric(),
  Lat = numeric(),
  Warmest_Month_SST = numeric(),
  SE = numeric(),
  Mask = numeric(),
  stringsAsFactors = FALSE
)

#loop through each grid cell
for (i in lon_indices) {
  for (j in lat_indices) {
    mask_value <- mask[i, j]
    if (is.na(mask_value) || mask_value == 0) next
    
    #extract the warmest month index and its SST value
    warmest_month_index <- mmm_month[i, j]
    if (is.na(warmest_month_index)) next  #skip if no data
    
    warmest_month_sst <- clim_monthly[i, j, warmest_month_index]
    if (is.na(warmest_month_sst)) next  #skip if no data
    
    #extract SST values for all months at this grid cell
    monthly_sst <- clim_monthly[i, j, ]
    
    #skip if no monthly data
    if (all(is.na(monthly_sst))) next
    
    #calculate standard error
    std_dev <- sd(monthly_sst, na.rm = TRUE)
    se <- std_dev / sqrt(length(monthly_sst))
    
    #append to results
    results <- rbind(results, data.frame(
      Lon = lon[i],
      Lat = lat[j],
      Warmest_Month_SST = warmest_month_sst,
      SE = se,
      Mask = mask_value
    ))
  }
}

#save results to a CSV file
write.csv(results, "warmest_month_sst_bounding_box.csv", row.names = FALSE)

#close the NetCDF file
nc_close(nc_data)
