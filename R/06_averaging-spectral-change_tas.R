## takes all tas GCM spectral exponent data and calculates average across 
path = "/Volumes/NIKKI/CMIP5-GCMs/" 

## create vector of file folders to put data into:
gcm_models <- c("01_CMCC-CMS_tas",
                "02_GFDL-CM3_tas",
                "03_GFDL-ESM2G_tas",
                "04_HadGEM2-ES_tas",
                "05_inmcm4_tas",
                "06_IPSL-CM5A-MR_tas",
                "07_MIROC-ESM-CHEM_tas",
                "08_MIROC5_tas",
                "09_MPI-ESM-LR_tas",
                "10_MPI-ESM-MR_tas",
                "11_MRI-CGCM3_tas")

folders <- paste(path, gcm_models, "/", sep = "")

types = c("l_estimate_PSD_low", "s_estimate_PSD_low", "l_estimate_PSD_high", "s_estimate_PSD_high",
          "l_estimate_AWC", "s_estimate_AWC")


## loop through different sensitivity layers 
i = 1
while (i < 6) {
  
  gcm = 1
  ## loop through gcms
  while(gcm <= 11) {

    if(gcm == 1) {
      big_boi = readRDS(paste(folders[gcm], gcm_models[gcm], "_spectral-change_multifrac.rds", sep = ""))[[i]]
    }
    else {
      big_boi = stack(big_boi,
                      readRDS(paste(folders[gcm], gcm_models[gcm], "_spectral-change_multifrac.rds", sep = ""))[[i]])
    }
    
    gcm = gcm + 1
  }
  
  ## save 
  saveRDS(big_boi, paste(path, "ALL/all-gcms_", types[[i]], ".rds", sep = ""))
  
  i = i + 1
}