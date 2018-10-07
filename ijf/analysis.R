# devtools::install_github("carlanetto/M4comp2018")
library(M4comp2018)
library(pbapply)
library(forecastHybrid)
library(xtable)
library(ggplot2)
library(reshape2)

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

# Clean the names of the models for the LaTeX tables
cleanNames <- function(x){
    x[x == "ARIMA"] <- "Arima"
    x[x == "bats"] <- "BATS"
    x[x == "thetam"] <- "Theta"
    return(x)
}

load("bestMods.RData")
Selection <- cleanNames(sapply(bestMods, function(x) class(x)[1]))
selectionDistribution <- table(Selection)

load("oracleMods.RData")
Oracle <- cleanNames(sapply(oracleMods, function(x) class(x)[1]))
oracleModDistribution <- table(Oracle)

# Confusion matrix of selected model vs oracle
cm <- confusionMatrix(data=factor(Selection), reference=factor(Oracle))
xtab <- xtable(cm$table, caption="Reference best model on test set vs model from selection procedure")
addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{3}{c}{Reference} \\\\\n",
"Selected & Arima & BATS & Theta  \\\\n")
print(xtable(xtab), add.to.row = addtorow, include.colnames = FALSE)

# Sensitivity, specificity, etc
xtable(cm$byClass)

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
rm(ensembleForecasts, selectionForecasts, oracleForecasts)

MASE <- data.frame(selection = selectionMASE, ensemble = ensembleMASE)# , oracle = oracleMASE
MASE <- melt(MASE)
MASE$MASE <- MASE$value
MASE$model <- MASE$variable
MASE$value <- MASE$variable <- NULL

####################################################################################################
# Distributions
####################################################################################################
png("distribution.png", width = 1200, height = 720)
ggplot(MASE, aes(x=log(MASE), fill = model)) + geom_density(alpha=0.5)
dev.off()
