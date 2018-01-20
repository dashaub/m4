library(data.table)
library(pbapply)
library(compiler)
library(Mcomp)
library(doMC)
library(Tcomp)


cleanM <- function(mObj){
    for(i in seq_along(mObj)){
        if(grepl("DEMOGR", mObj[[i]]$type, ignore.case = TRUE)){
            mObj[[i]]$type <- "DEMOGRAPHIC"
        } else if(grepl("INDUST", mObj[[i]]$type, ignore.case = TRUE)){
            mObj[[i]]$type <- "INDUSTRY"
        } else if(grepl("MACRO", mObj[[i]]$type, ignore.case = TRUE)){
            mObj[[i]]$type <- "MACRO"
        } else if(grepl("MICRO", mObj[[i]]$type, ignore.case = TRUE)){
            mObj[[i]]$type <- "MICRO"
        } else if(mObj[[i]]$type == "TOURISM"){
            mObj[[i]]$type <- "MACRO"
        } #else if(mObj[[i]]$period == "OTHER"){
            #mObj[[i]]$period <- "DAILY"
            #}
        }
    return(mObj)
    }

data(M3)
data(M1)
allData <- c(M1, M3, tourism)

createMObject <- function(x, type){
    xClean <- na.interp(x)
    tsLength <- length(x)
    # Ensure msts for hourlyl
    if(frequency(xClean) == 8766){
        xClean <- msts(xClean, seasonal.periods = c(24, 168, 8766), ts.frequency = 24)
        }
    # Ensure msts for daily
    if(frequency(xClean) == 365){
        xClean <- msts(xClean, seasonal.periods = c(7, 365.25), ts.frequency = 7)
        }
    horizon <- getHorizonFromFrequency(xClean)
    testSet <- subset(xClean, start = tsLength - horizon + 1)
    trainSet <- subset(xClean, end = tsLength - horizon)
    returnList <- list(x = trainSet, xx = testSet, h = horizon, n = length(trainSet), type = type)
    return(returnList)
    }

if(require(tscompdata)){
    nn3Clean <- lapply(nn3, FUN = function(x) createMObject(x, type = "MICRO"))
    nn5Clean <- lapply(nn5, FUN = function(x) createMObject(x, type = "MICRO"))
    gefcom2012_loadClean <- lapply(gefcom2012_load,
                                   FUN = function(x) createMObject(x, type = "MICRO"))
    gefcom2012_tempClean <- lapply(gefcom2012_temp,
                                   FUN = function(x) createMObject(x, type = "OTHER"))
    gefcom2012_wpClean <- lapply(gefcom2012_wp,
                                 FUN = function(x) createMObject(x, type = "MICRO"))
    allData <- c(allData, nn3Clean, nn5Clean,
                 gefcom2012_loadClean, gefcom2012_tempClean, gefcom2012_wpClean)
    
    }


cleaned <- cleanM(allData)
# Shuffle data so it distributes evenly for parallel feature extraction
set.seed(31415926)
cleaned <- rclean <- sample(cleaned, size = length(cleaned), replace = FALSE)
# Remove short series
shortSeries <- sapply(cleaned, FUN = function(x) length(x$x) <= 9)
cleaned <- cleaned[!shortSeries]
longSeries <- sapply(cleaned, FUN = function(x) length(x$x) > 5000)
cleaned <- cleaned[!longSeries]
save(cleaned, file = "cleaned.RData", compress = "xz", compression_level = 9)
