# devtools::install_github("carlanetto/M4comp2018")
library(M4comp2018)
library(pbapply)
library(forecastHybrid)

numCores <- 2

# Full ensemble forecasts
load("mods.RData")
ensembleForecasts <- pblapply(mods, function(x) forecast(x), cl = numCores)
rm(mods)
gc()

save(ensembleForecasts, file = "ensembleForecasts.RData")

