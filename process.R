library(data.table)
library(forecastHybrid)
library(tseries)
library(e1071)
library(uroot)
library(Mcomp)

# Load all the data files
files <- dir(pattern = "*.csv")


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
            }
        }
    return(mObj)
    }

data(M3)
data(M1)
allData <- c(M1, M3)

# All possible models
models <- c("a", "e", "f", "n", "s", "t")
expandedGrid <- expand.grid(rep(list(models), times = length(models)))
noDupes <- unique(apply(expandedGrid, 1, FUN = function(x) sort(unique(x))))
allModels <- noDupes[sapply(noDupes, FUN = function(x) length(x) >= 2)]
allModels <- sapply(allModels, FUN = function(x) paste0(x, collapse=""))

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
                    n_s_diff = arima_order[7]
                    )
    options(warn=0)
    return(df)
    }

cleaned <- cleanM(allData)

determineType <- function(df){
    # For each data file, determine if it is daily, hourly, weekly, etc
    return(substr(head(df$V1, 1), start = 0, 1))
    }

getHorizon <- function(x){
    # Determine the forecast horizon
    horizon <- c("H" = 48, "D" = 14, "W" = 13, "M" = 18, "Y" = 6)
    return(horizon[x])
    }

processFile <- function(x){
    # Load the data, train a model, produce forecasts, and write results
    dat <- fread(x, header = TRUE)
    type <- determineType(dat)
    dat$V1 <- NULL
    horizon <- getHorizon(type)
    dat = dat[, lapply(.SD, as.numeric)]
    return(dat)
    }


################################################################################
set.seed(314159265)
someSeries <- sample(cleaned, 50)
library(caret)
library(ranger)
# Prepare data for training
