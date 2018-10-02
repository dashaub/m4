library(forecastHybrid)
library(data.table)
library(pbapply)

currentSeries <- "Yearly"
numCores <- 4


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


getHorizon <- function(input){
  tab <- c("Hourly" = 48, "Daily" = 14, "Weekly" = 13, "Monthly" = 18,
           "Quarterly" = 8, "Yearly" = 6)
  return(as.numeric(tab[input]))
}

message("Processing ", currentSeries)
inputPath <- paste0("~/m4/Data/", currentSeries, "-train.csv")
dat <- fread(inputPath, header = TRUE, data.table = FALSE)
# Extract series names and remove from dataframe
seriesNames <- dat[, 1]
dat <- dat[, -1]
h <- getHorizon(currentSeries)

# Transform to list
dat <- apply(dat, MARGIN = 1, FUN = function(x) extractList(x, currentSeries))
names(dat) <- seriesNames


# Fit the ensemble model. This include both the model used in M4 forecasting as well as the
# individual component models that will be compared. These are the final fit models
models <- ifelse(currentSeries == "Monthly", "fs", "aft")
mods <- pblapply(dat, function(x) hybridModel(x, models = models, verbose = FALSE), cl = numCores)
save(mods, file = "mods.RData")


# Fit again on the subset of the data. This is the data that will be used along with the holdout
# set for model selection
subsetMods <- pblapply(dat, function(x) hybridModel(tail(x, -h), models = models, verbose = FALSE),
                       cl = numCores)
save(subsetMods, file = "subsetMods.RData")


# Chose a single "best" model based on the holdout set
chooseBestModel <- function(x){
  mod <- x$mod
  data <- x$data
  # Create a holdout set
  holdout <- as.numeric(tail(data, h))

  bestModel <- NULL
  bestMASE <- Inf
  for(model in c("auto.arima", "thetam", "tbats")){
    prediction <- forecast(mod[[model]], h=h)
    currentMASE <- accuracy(prediction, holdout)["Test set", "MASE"]
    if(currentMASE < bestMASE){
      bestMase <- currentMASE
      bestModel <- mod[[model]]
    }
  }
  return(bestModel)
}

# Combine the model and data into a single list so we can easily parallelize
packageModelsAndData <- function(models, data){
  output <- vector(mode = "list", length = length(models))
  for(i in seq_along(models)){
    output[[i]] <- list(model = models[[i]], data = data[[i]])
  }
  return(output)
}

modelsAndData <- packageModelsAndData(subsetMods, dat)
rm(subsetMods)
gc()

bestMods <- pblapply(modelsAndData, function(x) chooseBestModel(x), cl = numCores)
save(bestMods, file = "bestMods.RData")
