## modifying wget batch scripts to download gcms to compute canada computer
library(tidyverse)

dir = "data-raw/wget_scripts/wget_tas"

sh_files <- c("wget_01_CMCC-CESM_historical.sh",
              "wget_01_CMCC-CESM_rcp85.sh",
              "wget_03_CMCC-CMS_historical.sh",
              "wget_03_CMCC-CMS_rcp85.sh",
              "wget_04_MPI-ESM-LR_historical.sh",
              "wget_04_MPI-ESM-LR_rcp85.sh",
              "wget_05_MPI-ESM-MR_historical.sh",
              "wget_05_MPI-ESM-MR_rcp85.sh",
              "wget_06_GFDL-ESM2G_historical.sh",
              "wget_06_GFDL-ESM2G_rcp85.sh",
              "wget_07_GFDL-CM3_historical.sh",
              "wget_07_GFDL-CM3_rcp85.sh",
              "wget_08_GFDL-ESM2M_historical.sh",
              "wget_08_GFDL-ESM2M_rcp85.sh",
              "wget_09_HadGEM2-CC_historical.sh",
              "wget_09_HadGEM2-CC_rcp85.sh",
              "wget_10_HadGEM2-ES_historical.sh",
              "wget_10_HadGEM2-ES_rcp85.sh",
              "wget_11_HadGEM2-AO_historical.sh",
              "wget_11_HadGEM2-AO_rcp85.sh",
              "wget_12_IPSL-CM5A-LR_historical.sh",
              "wget_12_IPSL-CM5A-LR_rcp85.sh",
              "wget_13_IPSL-CM5B-LR_historical.sh",
              "wget_13_IPSL-CM5B-LR_rcp85.sh",
              "wget_14_MIROC5_historical.sh",
              "wget_14_MIROC5_rcp85.sh",
              "wget_15_MIROC-ESM-CHEM_historical.sh",
              "wget_15_MIROC-ESM-CHEM_rcp85.sh",
              "wget_16_MIROC-ESM_historical.sh",
              "wget_16_MIROC-ESM_rcp85.sh",
              "wget_17_inmcm4_historical.sh",
              "wget_17_inmcm4_rcp85.sh",
              "wget_18_CNRM-CM5_historical.sh",
              "wget_18_CNRM-CM5_rcp85.sh",
              "wget_19_MRI-CGCM3_historical.sh",
              "wget_19_MRI-CGCM3_rcp85.sh",
              "wget_20_MRI-ESM1_historical.sh",
              "wget_20_MRI-ESM1_rcp85.sh",
              "wget_21_IPSL-CM5A-MR_historical.sh",
              "wget_21_IPSL-CM5A-MR_rcp85.sh")

link_codes <- c(rep("http://aims3.llnl.gov/thredds/fileServer/", 4), #1,3
                rep("http://esgf.nci.org.au/thredds/fileServer/", 2), #4
                rep("http://aims3.llnl.gov/thredds/fileServer/", 2), #5
                rep("http://esgdata.gfdl.noaa.gov/thredds/fileServer/gfdl_dataroot/", 6), #6,7,8
                rep("http://esgf-data1.ceda.ac.uk/thredds/fileServer/", 4), #9,10
                rep("http://aims3.llnl.gov/thredds/fileServer",2), #11
                rep("http://vesg.ipsl.upmc.fr/thredds/fileServer", 4),#12,13
                rep("http://esgf-data1.diasjp.net/thredds/fileServer/", 6),#14,15,16
                rep("http://aims3.llnl.gov/thredds/fileServer/", 2), #17
                rep("http://esg1.umr-cnrm.fr/thredds/fileServer", 2), #18
                rep("http://esgf-data1.diasjp.net/thredds/fileServer/", 4),#19,20
                rep("http://vesg.ipsl.upmc.fr/thredds/fileServer", 2))#21

## read files in one by one
for (i in 1:length(sh_files)) {
  sh <- readLines(paste(dir, sh_files[i], sep=""))
  sh <- paste(sh, collapse = ' ')
  
  ## extract just the links by splitting
  split = strsplit(sh, regex("'", literal=TRUE), fixed=TRUE)[[1]]

  links = split[which(str_detect(split, link_codes[i]))]
  
  ## turn into wget commands:
  links = paste("wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! '", links, "'\n", sep="")
  links = paste(links, collapse = "\n")
  
  ## now add batch script stuff:
  sbatch = "#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  "
  
  new_sh <- factor(paste(sbatch, links, sep = "\n"))
  
  write.table(new_sh, file = paste("data-processed/wget_scripts/", sh_files[i], sep = ''),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


dir = "data-raw/wget_sst/"

sh_files <- c("wget_01_CMCC-CESM_historical.sh",
              "wget_01_CMCC-CESM_rcp85.sh",
              "wget_03_CMCC-CMS_historical.sh",
              "wget_03_CMCC-CMS_rcp85.sh",
              "wget_04_MPI-ESM-LR_historical.sh",
              "wget_04_MPI-ESM-LR_rcp85.sh",
              "wget_05_MPI-ESM-MR_historical.sh",
              "wget_05_MPI-ESM-MR_rcp85.sh",
              "wget_06_GFDL-ESM2G_historical.sh",
              "wget_06_GFDL-ESM2G_rcp85.sh",
              "wget_07_GFDL-CM3_historical.sh",
              "wget_07_GFDL-CM3_rcp85.sh",
              "wget_08_GFDL-ESM2M_historical.sh",
              "wget_08_GFDL-ESM2M_rcp85.sh",
              "wget_09_HadGEM2-CC_historical.sh",
              "wget_09_HadGEM2-CC_rcp85.sh",
              "wget_10_HadGEM2-ES_historical.sh",
              "wget_10_HadGEM2-ES_rcp85.sh",
              "wget_11_HadGEM2-AO_historical.sh",
              "wget_11_HadGEM2-AO_rcp85.sh",
              "wget_12_IPSL-CM5A-LR_historical.sh",
              "wget_12_IPSL-CM5A-LR_rcp85.sh",
              "wget_13_IPSL-CM5B-LR_historical.sh",
              "wget_13_IPSL-CM5B-LR_rcp85.sh",
              "wget_14_MIROC5_historical.sh",
              "wget_14_MIROC5_rcp85.sh",
              "wget_15_MIROC-ESM-CHEM_historical.sh",
              "wget_15_MIROC-ESM-CHEM_rcp85.sh",
              "wget_16_MIROC-ESM_historical.sh",
              "wget_16_MIROC-ESM_rcp85.sh",
              "wget_17_inmcm4_historical.sh",
              "wget_17_inmcm4_rcp85.sh",
              "wget_18_CNRM-CM5_historical.sh",
              "wget_18_CNRM-CM5_rcp85.sh",
              "wget_19_MRI-CGCM3_historical.sh",
              "wget_19_MRI-CGCM3_rcp85.sh",
              "wget_20_MRI-ESM1_historical.sh",
              "wget_20_MRI-ESM1_rcp85.sh",
              "wget_21_IPSL-CM5A-MR_historical.sh",
              "wget_21_IPSL-CM5A-MR_rcp85.sh")

link_codes <- c(rep("http://esgf-node.cmcc.it/thredds/fileServer/", 4), #1,3
                rep("http://esgf.nci.org.au/thredds/fileServer/", 2), #4
                rep("http://aims3.llnl.gov/thredds/fileServer/", 2), #5
                rep("http://esgdata.gfdl.noaa.gov/thredds/fileServer/gfdl_dataroot/", 6), #6,7,8
                rep("http://esgf-data1.ceda.ac.uk/thredds/fileServer/", 4), #9,10
                rep("http://aims3.llnl.gov/thredds/fileServer",2), #11
                rep("http://vesg.ipsl.upmc.fr/thredds/fileServer", 4),#12,13
                rep("http://esgf-data1.diasjp.net/thredds/fileServer/", 6),#14,15,16
                rep("http://aims3.llnl.gov/thredds/fileServer/", 2), #17
                rep("http://esg1.umr-cnrm.fr/thredds/fileServer", 2), #18
                rep("http://esgf-data1.diasjp.net/thredds/fileServer/", 4),#19,20
                rep("http://vesg.ipsl.upmc.fr/thredds/fileServer", 2))#21

## read files in one by one
for (i in 1:length(sh_files)) {
  sh <- readLines(paste(dir, sh_files[i], sep=""))
  sh <- paste(sh, collapse = ' ')
  
  ## extract just the links by splitting
  split = strsplit(sh, regex("'", literal=TRUE), fixed=TRUE)[[1]]
  
  links = split[which(str_detect(split, link_codes[i]))]
  
  ## turn into wget commands:
  links = paste("wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! '", links, "'\n", sep="")
  links = paste(links, collapse = "\n")
  
  ## now add batch script stuff:
  sbatch = "#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  "
  
  new_sh <- factor(paste(sbatch, links, sep = "\n"))
  
  write.table(new_sh, file = paste("data-processed/wget_scripts/wget_tas", sh_files[i], sep = ''),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}



