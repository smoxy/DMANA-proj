###################################LOAD#########################################
load(url("https://github.com/smoxy/DMANA-proj/blob/main/MCA/Dati_DimRed.RData?raw=true"))
source("https://raw.githubusercontent.com/MatteoFasulo/Rare-Earth/main/R/util/coreFunctions.R")
loadPackages(c('ggplot2'))
rm(mtcars)
rm(USArrests)
laws <- Data_lawsB[,-c(2,7,12,17,23,29,34,40,45:51)]
################################################################################
