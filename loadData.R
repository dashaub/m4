library(thief)
library(forecastHybrid)
library(pbapply)
library(data.table)
library(parallel)
numCores <- 1

inputs <- c("Hourly", "Daily", "Weekly", "Monthly", "Quarterly", "Yearly")
paths <- paste0("~/m4/Data/", inputs,  "-train.csv")


getHorizon <- function(input){
  tab <- c("Hourly" = 48, "Daily" = 14, "Weekly" = 13, "Monthly" = 18, "Quarterly" = 8, "Yearly" = 6)
  return(as.numeric(tab[input]))
}

combineForecasts <- function(forecastPI, forecastsPoint){
  results <- forecastPI
  for(ind in seq_along(results)){
    extracted <- lapply(forecastsPoint, function(x) x[[ind]]$mean)
    extracted <- extracted[!sapply(extracted, is.null)]
    combined <- ts(rowMeans(data.frame(extracted), na.rm = TRUE))
    tsp(combined) <- tsp(forecastPI[[ind]]$mean)
    results[[ind]]$mean <- combined
    }
    return(results)
}

# Extract the dataframe to a list of msts objects
extractList <- function(x, seriesName){
  tab <- c("Hourly" = 24, "Daily" = 7, "Weekly" = 52, "Monthly" = 12, "Quarterly" = 4, "Yearly" = 1)
  seriesFrequency <- tab[seriesName]
  tab <- list("Hourly" = c(24, 168), "Daily" = c(7, 365.25))
  multSeason <- tab[[seriesName]]
  cleaned <- as.numeric(x[!is.na(as.numeric(x))])
  series <- msts(cleaned, seasonal.periods = multSeason, ts.frequency = seriesFrequency)
  return(series)
}

writeResults <- function(forecastList, seriesName){
  components <- c("mean", "lower", "upper")
  for(component in components){
    baseDir <- paste0("~/m4/", component)
    dir.create(file.path(baseDir), showWarnings = FALSE)
    forecasts <- t(data.frame(lapply(forecastList, FUN = function(x) x[[component]])))
    rNames <- rownames(forecasts)
    filename <- paste0(baseDir, seriesName, "-", component, ".csv")
    write.table(x = forecasts, file = filename, quote = FALSE, sep = ",", col.names = FALSE)
  }
}

# Temporary for debug
# Completed hourly, yearly, daily (weekly need at least 2 periods, and quarterly fail, monthly succeed with n=200 and single core)
allData <- inputs[5]
currentSeries <- allData
for(currentSeries in allData){
  message("Processing", currentSeries)
  inputPath <- paste0("~/m4/Data/", currentSeries, "-train.csv")
  dat <- fread(inputPath, header = TRUE, data.table = FALSE)
  # Extract series names and remove from dataframe
  seriesNames <- dat[, 1]
  dat <- dat[, -1]
  seriesHorizon <- getHorizon(currentSeries)

  # Transform to list
  dat <- apply(dat, MARGIN = 1, FUN = function(x) extractList(x, currentSeries))
  names(dat) <- seriesNames
  set.seed(1234)
  #dat <- sample(dat, 359)
  gc()

  # Generate the base forecasts for prediction intervals
  forecasts <- pblapply(dat,
                        FUN = function(x) forecast(hybridModel(y = x, models = "aft",
                                                               verbose = FALSE),
                                                   h = seriesHorizon, level = 95,
                                                   PI.combination = "mean"),
                        cl = numCores)

  if(currentSeries != "Yearly"){
    # Create point forecasts from an ensemble
    arimaForecasts <- pblapply(X = dat,
                               FUN = function(x) thief(x, h = seriesHorizon, usemodel = "arima"),
                               cl = numCores)
    thetaForecasts <- pblapply(X = dat,
                               FUN = function(x) thief(x, h = seriesHorizon, usemodel = "theta"),
                               cl = numCores)
    # Combine the forecasts
    combinedForecasts <- combineForecasts(forecasts, list(arimaForecasts, thetaForecasts))
  }

  # Write results
  writeResults(forecasts, currentSeries)
}
