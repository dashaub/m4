library(forecastHybrid)

library(caret)
library(ranger)
library(Boruta)

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


# Prepare data for training

#~ cleaned <- sample(allData, 50)
#~ table(sapply(cleaned, FUN = function(x) x$type))
#~ table(sapply(cleaned, FUN = function(x) x$period))
#~ summary(sapply(cleaned, FUN = function(x) length(x$x)))

# Fit multiple models and return the MASE
fitModels <- function(x, models, lambda = FALSE){
    results <- sapply(models, FUN = function(model) fitModel(x, method = model, lambda = lambda))
    df <- data.frame(matrix(results, nrow = 1))
    names(df) <- names(results)
    return(df)
    }
fitModels <- cmpfun(fitModels, options = list(optimize = 3))

# Fit a single model and return the MASE
fitModel <- function(series, method, lambda = FALSE){
    # Set lambda
    if(lambda){
        lambda <- BoxCox.lambda(series$x)
    } else{
        lambda <- NULL
    }
    # Too short for slm
    if(grepl("s", method) && frequency(series$xx) == 1){
        NA
    } else if(length(series$x) / frequency(series$x) <= 2){
        # Too short
        NA
    } else if(method == "a"){
        mod <- auto.arima(series$x, lambda = lambda)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "e"){
        mod <- ets(series$x, lambda = lambda)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "f"){
        mod <- thetam(series$x)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "n"){
        mod <- nnetar(series$x, lambda = lambda, MaxNWts = 10000)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "s"){
        mod <- stlm(series$x, lambda = lambda)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    }else if(method == "t"){
        mod <- tbats(series$x, lambda = lambda)
        fc <- forecast(mod, h = length(series$xx))
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    } else if(method == "z"){
        fc <- snaive(series$x, h = length(series$xx), lambda = lambda)
        as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    }else{
        if(grepl("n", method)){
            maxwts <- list(MaxNWts = 10000)
        }else{
            maxwts <- NULL
            }
        mod <- hybridModel(series$x, model = method, lambda = lambda,
                           verbose = FALSE, n.args = maxwts)
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
                           FUN = function(x) fitModels(x = x, models = allModels, lambda = FALSE),
                           cl = cl))
mase_lambda <- rbindlist(pblapply(X = cleaned,
                                  FUN = function(x) fitModels(x = x, models = allModels, lambda = TRUE),
                                  cl = cl))
save(mase, file = "mase.RData")
save(mase_lambda, file = "mase_lambda.RData")

#library(doParallel)
#registerDoParallel(8)
#registerDoMC(8)


means <- colMeans(mase, na.rm = TRUE)
sorted <- apply(mase, 1, FUN = sort)

# Create labels
min_mase <- apply(mase, 1, FUN = function(x) min(x, na.rm = TRUE))
min_lambda_mase <- apply(mase_lambda, 1, FUN = function(x) min(x, na.rm = TRUE)) 
use_lambda <- factor(min_mase >= min_lambda_mase)
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
# 0.3353251
# 177?
lambdaMod <- train(x = dat, y = use_lambda,
                        method = "ranger", trControl = tc,tuneLength = numMod)

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
