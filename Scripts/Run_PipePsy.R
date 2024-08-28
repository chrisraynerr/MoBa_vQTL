remove(list = ls())
source("SetupLocalRprofile_PipePsy.R")
setwd("N:/durable/projects/crayner/PipePsy/")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... loading packages \n")

packages <- c("dplyr","data.table","tidyr","foreach","stringr")

installed_packages <-  packages %in% rownames(installed.packages(lib.loc=LIB))
if(any(installed_packages==F)){
  install.packages(packages[!installed_packages], dependencies = T, lib=LIB)
}

suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only=T, lib.loc=LIB))
)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... loading data dictionary \n")

datd   <- 
  fread("DataDictionary_ForPipe.csv") %>% 
  filter(measure=="CBCL")

SCALES <- unique(datd$measure)

# ######################################################################################
foreach(i=seq_along(SCALES)) %do% {
   cat("\n Running: ",paste0(SCALES[i],"\n"))
   system(paste0("Rscript --vanilla Data_Pipe_Scales.R --scale=",SCALES[i]))
  }

#######################################################################################
foreach(i=seq_along(SCALES)) %do% {
  system(paste0("Rscript Data_Pipe_xIRT.R --scale=",SCALES[i]))
}