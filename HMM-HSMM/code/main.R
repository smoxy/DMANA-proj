###################################LOAD#########################################
load(url("https://github.com/smoxy/DMANA-proj/blob/main/HMM-HSMM/code/Examples_L31.RData?raw=true"))
source("https://raw.githubusercontent.com/MatteoFasulo/Rare-Earth/main/R/util/coreFunctions.R")
loadPackages(c('ggplot2','gamlss','tidyverse','tidyquant','magrittr','tseries','MVN'))
rm(returns)
rm(pollution)
rm(stock.names)
################################################################################

histDist(test, family = "BCT", nbins=50) #?
################################################################################
StockReturns %<>%
  mutate(Date=as.Date(StockReturns$Date,format="%d/%m/%Y"))

StockReturns2 <- StockReturns %>%
                  select(SP500,NASDAQ,ESTX50,FTSE,Date) %>%
                  filter(Date > "2003/01/01" & Date < "2016/06/23")

dfStruct <- function(dataframe){
  dates <- dataframe$Date
  SP500 <- dataframe$SP500
  NASDAQ <- dataframe$NASDAQ
  FTSE <- dataframe$FTSE
  ESTX50 <- dataframe$ESTX50

  name <- rep(c("SP500","NASDAQ","FTSE","ESTX50"),times=c(nrow(dataframe),nrow(dataframe),nrow(dataframe),nrow(dataframe)))
  dates <- rep(dates, times=4)
  value <- rep(c(SP500,NASDAQ,FTSE,ESTX50))

  df <- data.frame("name"=name,
                       "date"=as.Date(dates,format="%d/%m/%Y"),
                       "value"=value)
  return(df)
}

Stocks <- dfStruct(StockReturns2)

Stocks %>%
  ggplot(., aes(x=date, y=value, group=name))+
    geom_line(size=.7, col="gray13") +
    scale_x_date(date_labels = "%Y", date_breaks="2 year") +
    facet_wrap(~ name, ncol = 1, scale = "free_y") +
    labs(title = "ESTX50, FTSE, NASDAQ & SP500 Chart",
         subtitle = "Multiple Stocks",
         y = "Return",
         x = "") +
    theme_tq()

Stocks %>%
  group_by(name) %>%
  summarise("Mean"=round(mean(value),3),
            "Std. dev"=sd(value),
            "Skewness"=skewness(value),
            "Kurtosis"=kurtosis(value, excess=T),
            "Jarque–Bera test (p-value)"=paste(round(as.double(jarque.bera.test(value)$statistic),1)," (",as.double(jarque.bera.test(value)$p.value),")",sep=""))

qqnorm(StockReturns2$NASDAQ)
qqnorm(StockReturns2$ESTX50)
mvn(StockReturns2[,-c(1,4,5)],mvnTest="mardia",multivariatePlot="qq")
