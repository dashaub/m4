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


currentSeries <- inputs[3]

writeResults <- function(x, seriesName){
  # Write point forecasts
  baseDir <- "~/m4/forecasts/"
  dir.create(file.path(baseDir), showWarnings = FALSE)
  pointForecasts <- t(data.frame(lapply(x, FUN = function(x) x$mean)))
  filename = paste0(baseDir, seriesName, "-point.csv")
  write.table(x = pointForecasts, file = filename, quote = FALSE, sep = ",", col.names = FALSE)
}

for(currentSeries in allData){
  inputPath <- paste0("~/m4/Data/", currentSeries, "-train.csv") 
  dat <- fread(inputPath, header = TRUE, data.table = FALSE)
  # Extract series names and remove from dataframe
  seriesNames <- dat[, 1]
  dat <- dat[, -1]
  seriesFrequency <- getFrequency(currentSeries)
  seriesHorizon <- getHorizon(currentSeries)

  # Extract the dataframe to a list
  extractList <- function(x){
    series <- as.numeric(x)
    series <- ts(series[!is.na(series)], f = seriesFrequency)
    return(series)
  }

  # Forecast function
  thiefForecast <- function(series){
    fc <- thief(y = series, h = seriesHorizon, usemodel = "arima")
    return(fc)
  }

  # Transform to list
  dat <- apply(dat, MARGIN = 1, FUN = function(x) extractList(x))
  names(dat) <- seriesNames
  gc()
  forecasts <- pblapply(X = dat, FUN = thiefForecast, cl = 6)

  # Write results
  writeResults(forecasts, currentSeries)
}
