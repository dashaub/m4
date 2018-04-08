library(thief)
library(pbapply)
library(data.table)

inputs <- c("Daily", "Hourly", "Monthly", "Yearly", "Weekly", "Quarterly")
paths <- paste0("~/m4/Data/", inputs,  "-train.csv")

getFrequency <- function(input){
  inputNames <- c("Daily", "Hourly", "Monthly", "Yearly", "Weekly", "Quarterly")
  mapping <- c(7, 24, 12, 1, 52, 4)
  names(mapping) <- inputNames
  return(as.numeric(mapping[input]))
}

getHorizon <- function(input){
  inputNames <- c("Daily", "Hourly", "Monthly", "Yearly", "Weekly", "Quarterly")
  mapping <- c(14, 48, 18, 6, 13, 8)
  names(mapping) <- inputNames
  return(as.numeric(mapping[input]))
}


currentSeries <- inputs[1]

for(currentSeries in allData){
  inputPath <- paste0("~/m4/Data/", currentSeries, "-train.csv") 
  #dat <- read.csv(inputPath, header = TRUE, quote = '"')
  dat <- fread(inputPath, header = TRUE, data.table = FALSE)
  seriesNames <- dat[, 1]
  dat <- dat[, -1]
  seriesFrequency <- getFrequency(currentSeries)
  seriesHorizon <- getHorizon(currentSeries)
  dat <- dat[, -1]
  #dat <- apply(X = dat, MARGIN = 2, FUN = function(x) as.numeric(x))
  fclist <- list()
  for(ind in 1:ncol(dat)){
    print(paste(ind, "of", ncol(dat)))
    rawSeries <- as.numeric(dat[, ind])
    series <- ts(rawSeries[!is.na(rawSeries)], f = seriesFrequency)
    fc <- thief(y = series, h = seriesHorizon, usemodel = "theta")
    fclist[[ind]] <- fc
  }
}


extractList <- function(x){
  series <- as.numeric(x)
  series <- ts(series[!is.na(series)], f = seriesFrequency)
  return(series)
  }

dat <- apply(dat, MARGIN = 1, FUN = function(x) extractList(x))
names(dat) <- seriesNames
gc()

thiefForecast <- function(series){
  fc <- thief(y = series, h = seriesHorizon, usemodel = "arima")
}
paste("Processing", currentSeries)
res <- pblapply(X = dat, FUN = thiefForecast, cl = 2)
