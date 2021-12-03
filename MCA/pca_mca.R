###################################LOAD#########################################
load(url("https://github.com/smoxy/DMANA-proj/blob/main/MCA/Dati_DimRed.RData?raw=true"))
source("https://raw.githubusercontent.com/MatteoFasulo/Rare-Earth/main/R/util/coreFunctions.R")
loadPackages(c('ggplot2','ISOcodes'))
rm(mtcars)
rm(USArrests)
laws <- Data_lawsB[,-c(2,7,12,17,23,29,34,40,45:51)]
laws <- as.data.frame(lapply(laws, factor))
################################################################################
whichNation<- function(ISOcode){
  ISOcode<-as.character(toupper(ISOcode))
  ISO<-ISOcodes::ISO_3166_1[,c(2,4)]
  nation<-ISO[which(ISO$Alpha_3==ISOcode),'Name']

  return(paste(ISOcode,":",nation))
}
