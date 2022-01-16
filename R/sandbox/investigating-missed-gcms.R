# investigating missed gcms
library(raster)

## ACCESS1.0: can use
data <- stack('data-raw/tas_day_ACCESS1-0_historical_r1i1p1_18500101-18741231.nc')

plot(data[[1]])
