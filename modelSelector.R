library(caret)
library(ranger)
library(data.table)
library(forecastHybrid)
library(tseries)
library(e1071)
library(uroot)
library(Mcomp)
library(doMC)
library(Tcomp)

# Load all the data files
files <- dir(pattern = "*.csv")

# All possible models
#models <- c("a", "e", "f", "n", "s", "t")
#expandedGrid <- expand.grid(rep(list(models), times = length(models)))
#noDupes <- unique(apply(expandedGrid, 1, FUN = function(x) sort(unique(x))))
#allModels <- noDupes[sapply(noDupes, FUN = function(x) length(x) >= 1)]
#allModels <- sapply(allModels, FUN = function(x) paste0(x, collapse=""))
#allModels <- c(allModels, "z")
allModels <- c("a", "e", "f", "n", "s", "t", "z")

tsFeatures <- function(x){
    options(warn=-1)
    adf <- adf.test(x)
    scaled <- scale(x)
    min_x = min(scaled)
    max_x = max(scaled)
    tsFreq <- frequency(x) > 1
    s_trend = ifelse(tsFreq,
                     as.numeric(abs(coef(tslm(x ~ trend + season, x))["trend"])),
                     as.numeric(abs(coef(tslm(x ~ trend, x))["trend"])))
    trend2 <- AIC(tslm(x ~ trend, x)) > AIC(tslm(x ~ trend + I(trend^2), x))
    pp <- pp.test(x)
    pp_explosive <- pp.test(x, alternative = "explosive")
    arima_order <- tryCatch({auto.arima(x, method = "CSS")$arma},
                            error = function(error_condition){
                                        auto.arima(x)$arma
                                      })
    df = data.frame(len = length(x),
                    ndiffs = ndiffs(x),
                    adf_statistic = as.numeric(adf$statistic),
                    adf_p = as.numeric(adf$p.value),
                    s_trend = s_trend,
                    trend2 = trend2,
                    var = var(x),
                    min = min_x,
                    max = max_x,
                    abs = max(abs(c(min_x, max_x))),
                    range = max(x) - min(x),
                    skew = abs(skewness(x)),
                    kurtosis = abs(kurtosis(x)),
                    bc_lambda = BoxCox.lambda(x),
                    pp = pp$p.value,
                    pp_explosive = pp_explosive$p.value,
                    ar = arima_order[1],
                    ma = arima_order[2],
                    sar = arima_order[3],
                    sma = arima_order[4],
                    n_diff = arima_order[6],
                    n_s_diff = arima_order[7],
                    num_periods = length(x) / frequency(x)
                    )
    options(warn=0)
    return(df)
    }

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
            }
        }
        if(mObj[[i]]$period == "OTHER"){
            mObj[[i]]$period <- "DAILY"
            }
    return(mObj)
    }

data(M3)
data(M1)
allData <- c(M1, M3, tourism)
cleaned <- cleanM(allData)

# Prepare data for training
set.seed(31415926)
#cleaned <- sample(allData, 50)
#table(sapply(cleaned, FUN = function(x) x$type))
#table(sapply(cleaned, FUN = function(x) x$period))
#summary(sapply(cleaned, FUN = function(x) length(x$x)))

features <- rbindlist(lapply(cleaned, FUN = function(x) tsFeatures(x$x)))
save(features, file = "features.RData")


registerDoMC(8)
set.seed(50)
#Evaluate accuracy of methods
mase <- matrix(NA, nrow = length(cleaned), ncol = length(allModels))
colnames(mase) <- allModels
count <- 1
numSeries <- length(cleaned)
for(series in cleaned){
    # Naive method that wastefully refits all models
    print(paste0(count, " of ", numSeries))
    res <- foreach(i = seq_along(allModels), .packages="forecastHybrid") %dopar%{
        method = allModels[i]
        # Nonseasonal
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
            mod <- nnetar(series$x)
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
            mod <- hybridModel(series$x, model = method, verbose = FALSE)
            fc <- forecast(mod, h = length(series$xx), PI = FALSE)
            as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
            }
        }
    mase[count, ] <- unlist(res)
    count <- count + 1
    }
#stopImplicitCluster()

orders <- apply(mase, 2, FUN = order)
means <- colMeans(mase, na.rm = TRUE)
#selectMods <- unique(c(allModels[order(means)][1:7], "aef", "at", "an", "aen", "aet"))
selectMods <- allModels
labels <- factor(selectMods[apply(mase[, selectMods], 1,  FUN = which.min)])
save(labels, file = "labels.RData")

dat <- features
dat$labels <- labels
dat$type <- as.character(sapply(cleaned, FUN = function(x) x$type))
set.seed(34)
seeds <- sample(1:(2*10^6), 2*10^6)
tc <- trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "random")
mod <- train(labels ~ ., data = dat, method = "ranger", trControl = tc, tuneLength = 3)
