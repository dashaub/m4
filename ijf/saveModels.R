library(forecastHybrid)
library(data.table)
library(pbapply)
if(!require(M4comp2018)){
  devtools::install_github("carlanetto/M4comp2018")
  library(M4comp2018)
}

currentSeries <- "Yearly"
numCores <- 6


####################################################################################################
# Prepare data
####################################################################################################


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


####################################################################################################
# Ensemble models
####################################################################################################


# Fit the ensemble model. This includes both the model used in M4 forecasting as well as the
# individual component models that will be compared. These are the final fit models
models <- ifelse(currentSeries == "Monthly", "fs", "aft")
mods <- pblapply(dat, function(x) hybridModel(x, models = models, verbose = FALSE), cl = numCores)
save(mods, file = "mods.RData")
gc()


####################################################################################################
# Helper functions for model selection
####################################################################################################


# Combine the M4 train and test data for yearly
combineTrainTestSet <- function(x){
  dat <- x[sapply(x, function(i) i$period == "Yearly")]
  res <- lapply(dat, function(i) c(i$x, i$xx))
  return(res)
}

# Combine the model and data into a single list so we can easily parallelize
packageModelsAndData <- function(models, data){
  output <- lapply(seq_along(models), function(x) list(model = models[[x]], data = data[[x]]))
  return(output)
}

# Chose a single "best" model based on the holdout set
chooseBestModel <- function(x){
  mod <- x$model
  data <- x$data
  # Create a holdout set
  holdout <- as.numeric(tail(data, h))

  bestModel <- NULL
  bestMASE <- Inf
  perfectModels <- list()
  for(model in c("auto.arima", "tbats", "thetam")){
    currentModel <- mod[[model]]
    prediction <- forecast(currentModel, h=h)
    currentMASE <- accuracy(prediction, holdout)["Test set", "MASE"]

    if(all(prediction$mean == holdout)){
      perfectModels[[model]] <- currentModel
      next
    }
    if(currentMASE < bestMASE){
      bestMASE <- currentMASE
      bestModel <- currentModel
    }
  }

  # If at least one model was a perfect fit, randomly select from those that fit perfectly
  if(length(perfectModels)){
    return(perfectModels[[sample(seq_along(perfectModels), 1)]])
    }
  return(bestModel)
}


####################################################################################################
# Oracle model
####################################################################################################

data(M4)
m4Data <- combineTrainTestSet(M4)
modelsAndData <- packageModelsAndData(mods, m4Data)
rm(mods, m4Data)
gc()

oracleMods <- pblapply(modelsAndData, function(x) chooseBestModel(x), cl = numCores)
save(oracleMods, file = "oracleMods.RData")
rm(oracleMods)
gc()


####################################################################################################
# Subset train data models
####################################################################################################


# Fit again on the subset of the data. This is the data that will be used along with the holdout
# set for model selection
subsetMods <- pblapply(dat, function(x) hybridModel(tail(x, -h), models = models, verbose = FALSE),
                       cl = numCores)
save(subsetMods, file = "subsetMods.RData")


####################################################################################################
# Model selection for best component model on holdout set
####################################################################################################


modelsAndData <- packageModelsAndData(subsetMods, dat)
rm(subsetMods)
gc()

bestMods <- pblapply(modelsAndData, function(x) chooseBestModel(x), cl = numCores)
save(bestMods, file = "bestMods.RData")
