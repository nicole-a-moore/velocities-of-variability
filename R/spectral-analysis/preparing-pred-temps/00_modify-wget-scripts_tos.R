## modifying wget batch scripts to download gcms to compute canada computer
library(tidyverse)

dir = "data-raw/wget_scripts/wget_tos/"

sh_files <- c("wget_01_CMCC-CMS_historical.sh",
              "wget_01_CMCC-CMS_rcp85.sh",
              "wget_02_GFDL-CM3_historical.sh",
              "wget_02_GFDL-CM3_rcp85.sh",
              "wget_03_GFDL-ESM2G_historical.sh",
              "wget_03_GFDL-ESM2G_rcp85.sh",
              "wget_04_HadGEM2-ES_historical.sh",
              "wget_04_HadGEM2-ES_rcp85.sh",
              "wget_05_INM-CM4_historical.sh",
              "wget_05_INM-CM4_rcp85.sh",
              "wget_06_IPSL-CM5A-MR_historical.sh",
              "wget_06_IPSL-CM5A-MR_rcp85.sh",
              "wget_07_MIROC-ESM-CHEM_historical.sh",
              "wget_07_MIROC-ESM-CHEM_rcp85.sh",
              "wget_08_MIROC5_historical.sh",
              "wget_08_MIROC5_rcp85.sh",
              "wget_09_MPI-ESM-LR_historical.sh",
              "wget_09_MPI-ESM-LR_rcp85.sh",
              "wget_10_MPI-ESM-MR_historical.sh",
              "wget_10_MPI-ESM-MR_rcp85.sh",
              "wget_11_MRI-CGCM3_historical.sh",
              "wget_11_MRI-CGCM3_rcp85.sh")

link_codes <- c(rep("http://aims3.llnl.gov/thredds/fileServer/", 2), #1
                rep("http://esgdata.gfdl.noaa.gov/thredds/fileServer/gfdl_dataroot/", 4), #2,3
                rep("http://esgf-data1.ceda.ac.uk/thredds/fileServer/", 2),#4
                rep("http://aims3.llnl.gov/thredds/fileServer/", 2), #5
                rep("http://vesg.ipsl.upmc.fr/thredds/fileServer/", 2), #6
                rep("http://esgf-data1.diasjp.net/thredds/fileServer/", 4), #7,8
                rep("http://esgf-data1.ceda.ac.uk/thredds/fileServer/", 2), # 9
                rep("http://aims3.llnl.gov/thredds/fileServer/", 2), # 10
                rep("http://esgf-data1.diasjp.net/thredds/fileServer/", 2) #11
                )


## read files in one by one
for (i in 1:length(sh_files)) {
  sh <- readLines(paste(dir, sh_files[i], sep=""))
  sh <- paste(sh, collapse = ' ')
  
  ## extract just the links by splitting
  split = strsplit(sh, regex("'", literal=TRUE), fixed=TRUE)[[1]]

  links = split[which(str_detect(split, link_codes[i]))]
  
  ## turn into wget commands:
  links = paste("wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password= '", links, "'\n", sep="")
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
  
  write.table(new_sh, file = paste("data-processed/wget_scripts/wget_tos/", sh_files[i], sep = ''),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

