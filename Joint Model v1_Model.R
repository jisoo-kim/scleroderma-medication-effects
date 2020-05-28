
##### Packages ##### 
Packages <- c("plyr","tidyr","ggplot2","data.table", "MCMCglmm")
lapply(Packages, library, character.only = TRUE)

##### Data Read-in ##### 
load("tmp_JK.Rdata")


##### Fit mcmcglmm to FVC, DLCO, MMF as outcome variables #####
MCMCglmm(cbind(), data = dat)
