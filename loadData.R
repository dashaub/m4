library(thief)
library(forecastHybrid)
library(pbapply)
library(data.table)
numCores <- 6

inputs <- c("Daily", "Hourly", "Monthly", "Yearly", "Weekly", "Quarterly")
paths <- paste0("~/m4/Data/", inputs,  "-train.csv")


getFrequency <- function(input){
  tab <- c("Hourly" = 24, "Daily" = 7, "Weekly" = 52, "Monthly" = 12, "Quarterly" = 4, "Yearly" = 1)
  return(as.numeric(tab[input]))
}

getHorizon <- function(input){
  tab <- c("Hourly" = 48, "Daily" = 14, "Weekly" = 13, "Monthly" = 18, "Quarterly" = 8, "Yearly" = 6)
  return(as.numeric(tab[input]))
}

combineForecasts <- function(x){
  results <- x[[1]]
  for(ind in seq_along(results)){
    extracted <- lapply(x, function(x) x[[ind]]$mean)
    combined <- ts(rowMeans(data.frame(extracted)))
    tsp(combined) <- tsp(extracted)[[1]]
    results[[ind]]$mean <- combined
    }
    return(results)
}

# Extract the dataframe to a list
extractList <- function(x, seriesFrequency){
  series <- as.numeric(x)
  series <- ts(series[!is.na(series)], f = seriesFrequency)
  return(series)
}

writeResults <- function(x, seriesName){
  # Write point forecasts
  baseDir <- "~/m4/forecasts/"
  dir.create(file.path(baseDir), showWarnings = FALSE)
  pointForecasts <- t(data.frame(lapply(x, FUN = function(x) x$mean)))
  rNames <- rownames(pointForecasts)
  filename = paste0(baseDir, seriesName, "-point.csv")
  write.table(x = pointForecasts, file = filename, quote = FALSE, sep = ",", col.names = FALSE)

  # Write upper forecasts
  baseDir <- "~/m4/upper/"
  dir.create(file.path(baseDir), showWarnings = FALSE)
  upperForecasts <- t(data.frame(lapply(x, FUN = function(x) x$upper)))
  rownames(upperForecasts) <- rNames
  filename = paste0(baseDir, seriesName, "-upper.csv")
  write.table(x = upperForecasts, file = filename, quote = FALSE, sep = ",", col.names = FALSE)

  # Write lower forecasts
  baseDir <- "~/m4/lower/"
  dir.create(file.path(baseDir), showWarnings = FALSE)
  lowerForecasts <- t(data.frame(lapply(x, FUN = function(x) x$lower)))
  rownames(lowerForecasts) <- rNames
  filename = paste0(baseDir, seriesName, "-lower.csv")
  write.table(x = lowerForecasts, file = filename, quote = FALSE, sep = ",", col.names = FALSE)
}

# Temporary for debug
# Completed hourly, yearly, testing quarterly (weekly needs at least 2 periods)
allData <- inputs[6]
currentSeries <- allData
for(currentSeries in allData){
  print(paste("Processing", currentSeries))
  inputPath <- paste0("~/m4/Data/", currentSeries, "-train.csv") 
  dat <- fread(inputPath, header = TRUE, data.table = FALSE)
  # Extract series names and remove from dataframe
  seriesNames <- dat[, 1]
  dat <- dat[, -1]
  seriesFrequency <- getFrequency(currentSeries)
  seriesHorizon <- getHorizon(currentSeries)

  # Transform to list
  dat <- apply(dat, MARGIN = 1, FUN = function(x) extractList(x, seriesFrequency))
  names(dat) <- seriesNames
  gc()

  # Generate the base forecasts for prediction intervals
  forecasts <- pblapply(dat,
                        FUN = function(x) forecast(hybridModel(y = x, models = "aft",
                                                               verbose = FALSE),
                                                   h = seriesHorizon, level = 95,
                                                   PI.combination = "mean"),
                        cl = numCores)

  if(currentSeries != "Yearly"){
    # Forecast for prediction intervals
    forecastsPI <- pblapply(dat,
                            FUN = function(x) forecast(hybridModel(y = x, models = "aft",
                                                                   verbose = FALSE),
                                                       h = seriesHorizon, level = 95),
                            cl = numCores)
    # Two forecasts to ensemble
    arimaForecasts <- pblapply(X = dat,
                               FUN = function(x) thief(x, h = seriesHorizon, usemodel = "arima"),
                               cl = numCores)
    thetaForecasts <- pblapply(X = dat,
                               FUN = function(x) thief(x, h = seriesHorizon, usemodel = "theta"),
                               cl = numCores)
    # Combine the forecasts
    forecasts <- combineForecasts(list(arimaForecasts, thetaForecasts))
  }

  # Write results
  writeResults(forecasts, currentSeries)
}
