# devtools::install_github("carlanetto/M4comp2018")
library(M4comp2018)
library(pbapply)
library(forecastHybrid)
library(xtable)
library(ggplot2)

numCores <- 2

####################################################################################################
# Forecasts
####################################################################################################


# Full ensemble forecasts
load("mods.RData")
ensembleForecasts <- pblapply(mods, function(x) forecast(x), cl = numCores)
save(ensembleForecßasts, file = "ensembleForecasts.RData")

# Forecasts from the best model according to the model selection procedure
load("bestMods.RData")
selectionForecasts <- pblapply(bestMods, function(x) forecast(x), cl = numCores)
save(selectionForecasts, file = "selectionForecasts.RData")

# Forecasts from the best model given by on Oracle on the test set
load("oracleMods.RData")
oracleForecasts <- pblapply(oracleMods, function(x) forecast(x), cl = numCores)
save(oracleForecasts, file = "oracleForecasts.RData")


####################################################################################################
# Tables
####################################################################################################


# LaTeX table of distribution of selected best individual models
load("bestMods.RData")
bestModDistribution <- table(sapply(bestMods, function(x) class(x)[1]))

# LaTeX table of distribution of selected best individual models from oracle
load("oracleMods.RData")
oracleModDistribution <- table(sapply(oracleMods, function(x) class(x)[1]))


####################################################################################################
# MASE calculation
####################################################################################################


rm(oracleMods, bestMods, mods)
gc()

data(M4)
# Calculuate MASE for a specific model
calculateMASE <- function(fc){
  dat <- M4[sapply(M4, function(i) i$period == "Yearly")]
  lambda <- function(i) accuracy(fc[[i]], as.numeric(dat[[i]]$xx))["Test set", "MASE"]
  res <- pblapply(seq_along(fc), lambda, cl = numCores)
  res <- unlist(res)
  return(res)
}

ensembleMASE <- calculateMASE(ensembleForecasts)
selectionMASE <- calculateMASE(selectionForecasts)
oracleMASE <- calculateMASE(oracleForecasts)
save(ensembleMASE, oracleMASE, selectionMASE, file="MASE.RData")

####################################################################################################
# Distributions
####################################################################################################

