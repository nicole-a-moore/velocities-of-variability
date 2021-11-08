### seeing what spatial resolutions we have 
library(tidyverse)

# data from https://portal.enes.org/data/enes-model-data/cmip5/resolution
# If two values are given for the latitude resolution of the ocean grid, 
# resolution is not constant. The first value is that for the equator, the second 
# for the poles (maximum for the two poles if different)

sums <- read.csv("data-raw/GCM_summaries.csv")

head(sums)

gcm_models <- c("CMCC-CESM", "CMCC-CM", 'CMCC-CMS', 'MPI-ESM-LR', 'MPI-ESM-MR',
                "GFDL-ESM2G", 'GFDL-CM3', 'GFDL-ESM2M', "HadGEM2-CC", 'HadGEM2-ES',
                "HadGEM2-AO", 'IPSL-CM5A-LR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM-CHEM',
                'MIROC-ESM', "INM-CM4", 'CNRM-CM5', "MRI-CGCM3", 'MRI-ESM1',
                'IPSL-CM5A-MR')

length(which(gcm_models %in% unique(sums$Model)))

mods <- filter(sums, Model %in% gcm_models)

# min resolution:
mods$Model[which(mods$Lon_land == min(mods$Lon_land))] # CMCC-CM is lowest spatial resolution
min(mods$Lon_land) # 0.7484 deg lon, 0.75 deg lat
