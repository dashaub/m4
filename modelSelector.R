library(caret)
library(ranger)
library(Boruta)
library(data.table)
library(forecastHybrid)
library(pbapply)
library(compiler)

library(Mcomp)
library(doMC)
library(Tcomp)


# Load all the data files
files <- dir(pattern = "*.csv")

# All possible models
#~ models <- c("a", "e", "f", "n", "s", "t")
#~ expandedGrid <- expand.grid(rep(list(models), times = length(models)))
#~ noDupes <- unique(apply(expandedGrid, 1, FUN = function(x) sort(unique(x))))
#~ allModels <- noDupes[sapply(noDupes, FUN = function(x) length(x) >= 1)]
#~ allModels <- sapply(allModels, FUN = function(x) paste0(x, collapse=""))
#~ allModels <- c(allModels, "z")
allModels <- c("a", "e", "f", "n", "s", "t", "z")


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
    gefcom2012_loadClean <- lapply(gefcom2012_load, FUN = function(x) createMObject(x, type = "MICRO"))
    gefcom2012_tempClean <- lapply(gefcom2012_temp, FUN = function(x) createMObject(x, type = "OTHER"))
    gefcom2012_wpClean <- lapply(gefcom2012_wp, FUN = function(x) createMObject(x, type = "MICRO"))
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

# Prepare data for training

#~ cleaned <- sample(allData, 50)
#~ table(sapply(cleaned, FUN = function(x) x$type))
#~ table(sapply(cleaned, FUN = function(x) x$period))
#~ summary(sapply(cleaned, FUN = function(x) length(x$x)))

# Fit multiple models and return the MASE
fitModels <- function(x, models){
    results <- sapply(models, FUN = function(model) fitModel(x, method=model))
    df <- data.frame(matrix(results, nrow = 1))
    names(df) <- names(results)
    return(df)
    }
fitModels <- cmpfun(fitModels, options = list(optimize = 3))

# Fit a single model and return the MASE
fitModel <- function(series, method){
    if(grepl("s", method) && frequency(series$xx) == 1){
        NA
    } else if(length(series$x) / frequency(series$x) <= 2){
        # Too short
        NA
    } else if(method == "a"){
        mod <- auto.arima(series$x)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "e"){
        mod <- ets(series$x)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "f"){
        mod <- thetam(series$x)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "n"){
        mod <- nnetar(series$x, MaxNWts = 10000)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "s"){
        mod <- stlm(series$x)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    }else if(method == "t"){
        mod <- tbats(series$x)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "z"){
        fc <- snaive(series$x, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    }else{
        if(grepl("n", method)){
            maxwts <- list(MaxNWts = 10000)
        }else{
            maxwts <- NULL
            }
        mod <- hybridModel(series$x, model = method, verbose = FALSE, n.args = maxwts)
        fc <- forecast(mod, h = length(series$xx), PI = FALSE)
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
        }
    }
fitModel <- cmpfun(fitModel, options = list(optimize = 3))

# Serial features for debugging
res <- list()
for(i in seq_along(cleaned)){
    print(i)
    x <- cleaned[[i]]$x
    res[[i]] <- tsFeatures(x)
}


# Build features and labels
cl <- makeForkCluster(8)
#features <- rbindlist(parLapplyLB(cl = cl, X = cleaned, fun = function(x) tsFeatures(x$x)))
features <- rbindlist(pblapply(X = cleaned, FUN = function(x) tsFeatures(x$x), cl = cl))
save(features, file = "features.RData", compress = "xz", compression_level = 9)

#res <- parLapplyLB(cl = cl, X = cleaned, fun = function(x) fitModels(x = x, models = allModels))
mase <- rbindlist(pblapply(X = cleaned,
                           FUN = function(x) fitModels(x = x, models = allModels),
                           cl = cl))
save(mase, file = "mase.RData")


#library(doParallel)
#registerDoParallel(8)
#registerDoMC(8)


means <- colMeans(mase, na.rm = TRUE)
sorted <- apply(mase, 1, FUN = sort)

# Create labels
labels_first <- factor(sapply(sorted, FUN = function(x) names(x[1])))
labels_second <- factor(sapply(sorted, FUN = function(x) names(x[2])))
labels_third <- factor(sapply(sorted, FUN = function(x) names(x[3])))
labels_second_worst <- factor(sapply(sorted, FUN = function(x) names(tail(x, 2)[1])))
labels_worst <- factor(sapply(sorted, FUN = function(x) names(tail(x, 1))))
save(labels_first, file = "labels_first.RData")
save(labels_second, file = "labels_second.RData")
save(labels_third, file = "labels_third.RData")
save(labels_second_worst, file = "labels_second_worst.RData")
save(labels_worst, file = "labels_worst.RData")


# Feature selection
dat <- features
dat$type <- as.character(sapply(cleaned, FUN = function(x) x$type))
maxRuns <- 100
set.seed(34)
bFirst <- Boruta(x = dat, y = labels_first, maxRuns = maxRuns, doTrace = 1)
save(bFirst, file = "b.RData", compress = "xz", compression_level = 9)
set.seed(34)
bSecond <- Boruta(x = dat, y = labels_second, maxRuns = maxRuns, doTrace = 1)
save(bSecond, file = "b.RData", compress = "xz", compression_level = 9)
set.seed(34)
bThird <- Boruta(x = dat, y = labels_third, maxRuns = maxRuns, doTrace = 1)
save(bThird, file = "b.RData", compress = "xz", compression_level = 9)
set.seed(34)
bWorst <- Boruta(x = dat, y = labels_worst, maxRuns = maxRuns, doTrace = 1)
save(bWorst, file = "b.RData", compress = "xz", compression_level = 9)
set.seed(34)
bSecondWorst <- Boruta(x = dat, y = labels_second_worst, maxRuns = maxRuns, doTrace = 1)
save(bSecondWorst, file = "b.RData", compress = "xz", compression_level = 9)

# Train models
#registerDoMC(1)
tc <- trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "random")
numMod <- 150
seed <- 50
set.seed(seed)
rangerModFirst <- train(x = dat, y = labels_first,
                        method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModFirst, file = "rangerModFirst.RData",
     compress = "xz", compression_level = 9)
set.seed(seed)
rangerModSecond <- train(x = dat, y = labels_second,
                         method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModSecond, file = "rangerModSecond.RData",
     compress = "xz", compression_level = 9)
set.seed(seed)
rangerModThird <- train(x = dat, y = labels_third,
                        method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModThird, file = "rangerModThird.RData",
     compress = "xz", compression_level = 9)
set.seed(seed)
rangerModSecondWorst <- train(x = dat, y = labels_second_worst,
                              method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModSecondWorst, file = "rangerModSecondWorst.RData",
     compress = "xz", compression_level = 9)
set.seed(seed)
rangerModWorst <- train(x = dat, y = labels_worst,
                        method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModWorst, file = "rangerModWorst.RData",
     compress = "xz", compression_level = 9)
