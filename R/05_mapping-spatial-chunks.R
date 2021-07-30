## visualizing the output:

path = "/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/"

filepath = paste(path, "sp_files.rds", sep = "")
names = readRDS(filepath)
names = str_replace_all(names, 'CMIP5-GCMs', '/Volumes/SundayLab/CMIP5-GCMs')

l_filenames <- str_replace_all(names, "spatial_temps", 'l-detrended')
s_filenames <-  str_replace_all(names, "spatial_temps", 's-detrended')

lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

makeMaps <- function(date) {

  lat_index = 1
  lon_index = 1
  count = 1
  while(count < length(l_filenames)+1) {
    
    num = which(paste(year, date, sep = '.') == dates)
    seq = seq(from = num, to = 83950, by = 365)
    
    ## retrieve spatial chunk from nc file 
    tas <- stack(names[count])[[seq]]
    l_detrended_tas <- stack(l_filenames[count])[[seq]]
    s_detrended_tas <- stack(s_filenames[count])[[seq]]
    
    tas[is.na(tas)] <- -999
    l_detrended_tas[is.na(l_detrended_tas)] <- -999
    s_detrended_tas[is.na(s_detrended_tas)] <- -999
    
    tas <- data.frame(rasterToPoints(tas))
    temp = tas$y
    tas$y = tas$x
    tas$x = temp
    
    l_detrended_tas <- data.frame(rasterToPoints(l_detrended_tas))
    temp = l_detrended_tas$y
    l_detrended_tas$y = l_detrended_tas$x
    l_detrended_tas$x = temp
    
    s_detrended_tas <- data.frame(rasterToPoints(s_detrended_tas))
    temp = s_detrended_tas$y
    s_detrended_tas$y = s_detrended_tas$x
    s_detrended_tas$x = temp
    
    key = data.frame("type" = 'lat', 'val' = unique(tas$y))
    key = rbind(key, data.frame("type" = 'lon', 'val' = sort(unique(tas$x))))
    
    l_key = data.frame("type" = 'lat', 'val' = unique(l_detrended_tas$y))
    l_key = rbind(l_key, data.frame("type" = 'lon', 'val' = sort(unique(l_detrended_tas$x))))
    
    s_key = data.frame("type" = 'lat', 'val' = unique(s_detrended_tas$y))
    s_key = rbind(s_key, data.frame("type" = 'lon', 'val' = sort(unique(s_detrended_tas$x))))
    
    ## get lat and lon of chunk:
    lt = lat[lat_index:(lat_index+59)]
    ln = lon[lon_index:(lon_index+59)]
    
    key$new = l_key$new = s_key$new = append(lt, ln)
    
    tas = left_join(tas, filter(key, type == 'lon'), by = c("x" = "val")) %>%
      mutate(x = new) %>% select(-new, -type)
    tas = left_join(tas, filter(key, type == 'lat'), by = c("y" = "val")) %>%
      mutate(y = new) %>% select(-new, -type) %>% rename('x' = y, 'y' = x)
    
    l_detrended_tas = left_join(l_detrended_tas, filter(key, type == 'lon'), by = c("x" = "val")) %>%
      mutate(x = new) %>% select(-new, -type)
    l_detrended_tas = left_join(l_detrended_tas, filter(key, type == 'lat'), by = c("y" = "val")) %>%
      mutate(y = new) %>% select(-new, -type) %>% rename('x' = y, 'y' = x)
    
    s_detrended_tas = left_join(s_detrended_tas, filter(key, type == 'lon'), by = c("x" = "val")) %>%
      mutate(x = new) %>% select(-new, -type)
    s_detrended_tas = left_join(s_detrended_tas, filter(key, type == 'lat'), by = c("y" = "val")) %>%
      mutate(y = new) %>% select(-new, -type) %>% rename('x' = y, 'y' = x)
    
    tas <- rasterFromXYZ(tas)
    l_detrended_tas <- rasterFromXYZ(l_detrended_tas)
    s_detrended_tas <- rasterFromXYZ(s_detrended_tas)
    
    tas[tas == -999] = NA
    l_detrended_tas[l_detrended_tas == -999] = NA
    s_detrended_tas[s_detrended_tas == -999] = NA
    
    ## merge together
    if (count == 1) {
      mosaic <- tas
      l_mosaic  <- l_detrended_tas
      s_mosaic  <- s_detrended_tas
    }
    else {
      mosaic <- mosaic(mosaic, tas, fun = mean)
      l_mosaic  <- mosaic(l_mosaic, l_detrended_tas, fun = mean)
      s_mosaic  <- mosaic(s_mosaic, s_detrended_tas, fun = mean)
    }
    
    if (count == 6) {
      lat_index <- lat_index + 60
      lon_index <- 1
    }
    else if (count == 12) {
      lat_index <- lat_index + 60
      lon_index <- 1
    }
    else {
      lon_index <- lon_index + 60
    }
    
    count = count + 1
  }
  
  return(list(mosaic, l_mosaic, s_mosaic))
}

test = makeMaps("05.28")

writeRaster(mosaic, "data-processed/mosaic.grd")
writeRaster(s_mosaic, "data-processed/s_mosaic.grd")
writeRaster(l_mosaic, "data-processed/l_mosaic.grd")

