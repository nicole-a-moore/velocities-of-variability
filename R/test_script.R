### test R script

vec <- c("first arg", "second arg", "third arg")

command_args <- commandArgs(trailingOnly = TRUE)

if (length(command_args)==0) {
  stop("At least one argument must be supplied!", call.=FALSE)
} else if (length(command_args)==1) {
  chosen_one = command_args[1]
}

print(command_args)

chosen_one <- as.numeric(chosen_one)

## try to get command line args
print("Hello world!", stdout())
print(paste("The chosen one is the ", as.character(chosen_one), " :-)", sep = ""), stdout())


print(class(chosen_one))

df <- data.frame(chosen_one)

write.csv(df, "CMIP5-GCMs/chosen_one.csv")





### cheching output:
path = "/Volumes/SundayLab/CMIP5-GCMs/" 
gcm_models <- c("01_CMCC-CESM", "02_CMCC-CM", '03_CMCC-CMS', '04_MPI-ESM-LR', '05_MPI-ESM-MR',
                "06_GFDL-ESM2G", '07_GFDL-CM3', '08_GFDL-ESM2M', '09_HadGEM2-CC', '10_HadGEM2-ES',
                "11_HadGEM2-AO", '12_IPSL-CM5A-LR', '13_IPSL-CM5B-LR', '14_MIROC5', '15_MIROC5-ESM-CHEM',
                '16_MIROC5-ESM', "17_inmcm4", '18_CNRM-CM5', "19_MRI-CGCM3", '20_MRI-ESM1',
                '21_IPSL-CM5A-MR')
folders <- paste(path, gcm_models, "/", sep = "")


file2 <- "resampled_tas_day_CMCC-CM_historical_r1i1p1_18700101-18701231.nc"
file3 <- "resampled_tas_day_CMCC-CMS_historical_r1i1p1_19000101-19091231.nc"
file4 <- "resampled_tas_day_MPI-ESM-LR_rcp85_r1i1p1_20800101-20891231.nc"
file5 <- "resampled_tas_day_MPI-ESM-MR_rcp85_r1i1p1_20900101-21001231.nc"

file = paste(folders[2], file2, sep = "")
gcm = stack(file)
plot(gcm[[1]])
gcm

file = paste(folders[3], file3, sep = "")
gcm = stack(file)
plot(gcm[[100]])
gcm

file = paste(folders[4], file4, sep = "")
gcm = stack(file)
plot(gcm[[10]])
gcm

file = paste(folders[5], file4, sep = "")
gcm = stack(file)
plot(gcm[[15]])
gcm







## look at files 
library(raster)
files <- c("l-detrended_lon-120-180_lat--30--90.nc", "s-detrended_lon-120-180_lat--30--90.nc", 
           "spatial_temps_lon-120-180_lat--30--90.nc")

path = "/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/"

filenames = paste(path, files, sep = "")


s <- stack(filenames[2])

plot(s[[1]])
