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
        if(grepl("DEMOGR", mObj[[i]]$type)){
            mObj[[i]]$type <- "DEMOGRAPHIC"
        } else if(grepl("INDUST", mObj[[i]]$type)){
            mObj[[i]]$type <- "INDUSTRY"
        } else if(grepl("MACRO", mObj[[i]]$type)){
            mObj[[i]]$type <- "MACRO"
        } else if(grepl("MICRO", mObj[[i]]$type)){
            mObj[[i]]$type <- "MICRO"
            }
        }
    return(mObj)
    }

data(M3)
data(M1)
allData <- c(M1, M3)

tsFeatures <- function(x){
    options(warn=-1)
    adf <- adf.test(x)
    scaled <- scale(x)
    min_x = min(scaled)
    max_x = max(scaled)
    s_trend = as.numeric(abs(coef(tslm(x ~ trend + season, x))["trend"]))
    tsFreq <- frequency(x)
    ch <- 0
    # Strange assignment behavior from ch.test does not work in ifelse()
    if(tsFreq){
        ch <- ch.test(x)
        }
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
                    var = var(scaled),
                    min = min_x,
                    max = max_x,
                    abs = max(abs(c(min_x, max_x))),
                    range = max(x) - min(x),
                    skew = skewness(x),
                    kurtosis = kurtosis(x),
                    bc_lambda = BoxCox.lambda(x),
                    ch_max = ifelse(tsFreq, max(ch$pvalues), 0),
                    ch_mean = ifelse(tsFreq, mean(ch$pvalues), 0),
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
