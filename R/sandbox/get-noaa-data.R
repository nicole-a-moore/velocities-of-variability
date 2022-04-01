## script for accessing data from NOAA without working directly with .nc files
library(rerddap)

## begin by finding out the Dataset ID of the data you want to access:
## https://coastwatch.pfeg.noaa.gov/erddap/index.html

## load the info for your dataset using your Dataset ID
info <- info("ncdcOisst21Agg_LonPM180")

## try accessing the data:
time_series <- griddap(info,
                       time = c("1981-09-01", "1991-01-01"), ## choose the times you want data between
                       latitude = c(40.125, 40.125), ## choose the maximum and minimum lat and lon you want data for
                       longitude = c(-160.125, -160.125))

## note: don't expect this function to finish running quickly if you are asking it for a large chunk of data -- it's normal for it to take a while
temps <- data.frame(temps = time_series$data$sst, time  = 1:length(time_series$data$sst))
ggplot(temps, aes(x = time, y = temps)) + geom_line() + theme_bw()



## troubleshooting:
## if you get the error message "Error in R_nc4_open: NetCDF: Unknown file format" when running griddap(), there was an issue accessing the data from ERDDAP using info() function
## to try getting the data again, you will need to run the function cache_delete_all() before re-running the info() and gridapp() functions





# Addendum: Plotting Sea Surface Temperature
#--------------------------------------------------

library(ggplot2)
library(dplyr)

# pull data from area
temps_spatial <- griddap(info,
                         time = c("2020-04-04", "2020-04-04"), ## choose the times you want data between
                         latitude = c(23, 50), ## choose the maximum and minimum lat and lon you want data for
                         longitude = c(-88, -51))


# plot
temps_spatial$data %>%
  ggplot(aes(x=lon,y=lat,fill=sst)) +
  geom_raster() +
  scale_fill_gradient(na.value = "white") +
  theme_bw() + 
  labs(x="Longitude",y="Latitude",fill = "Temp (°C)")