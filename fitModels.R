library(forecastHybrid)
library(thief)
library(prophet)
library(compiler)

# Fit multiple models and return the MASE
fitModels <- function(x, models, lambda = FALSE){
    results <- sapply(models, FUN = function(model) fitModel(x, method = model, lambda = lambda))
    df <- data.frame(matrix(results, nrow = 1))
    names(df) <- names(results)
    return(df)
    }
fitModels <- cmpfun(fitModels, options = list(optimize = 3))

getAccuracy <- function(fc, series){
    return(as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"]))
    }
getAccuracy <- cmpfun(getAccuracy, options = list(optimize = 3))

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
        return(NA)
    } else if(length(series$x) / frequency(series$x) <= 2){
        # Too short
        return(NA)
    } else if(method == "a"){
        mod <- auto.arima(series$x, lambda = lambda)
        fc <- forecast(mod, h = length(series$xx))
    } else if(method == "e"){
        mod <- ets(series$x, lambda = lambda)
        fc <- forecast(mod, h = length(series$xx))
    } else if(method == "f"){
        mod <- thetam(series$x)
        fc <- forecast(mod, h = length(series$xx))
    } else if(method == "n"){
        mod <- nnetar(series$x, lambda = lambda, MaxNWts = 10000)
        fc <- forecast(mod, h = length(series$xx))
    } else if(method == "s"){
        mod <- stlm(series$x, lambda = lambda)
        fc <- forecast(mod, h = length(series$xx))
    }else if(method == "t"){
        mod <- tbats(series$x, lambda = lambda)
        fc <- forecast(mod, h = length(series$xx))
    } else if(method == "z"){
        fc <- snaive(series$x, h = length(series$xx), lambda = lambda)
        #as.numeric(accuracy(fc, x = series$xx)["Test set", "MASE"])
    }else{
        if(grepl("n", method)){
            maxwts <- list(MaxNWts = 10000)
        }else{
            maxwts <- NULL
            }
        mod <- hybridModel(series$x, model = method, lambda = lambda,
                           verbose = FALSE, n.args = maxwts)
        fc <- forecast(mod, h = length(series$xx), PI = FALSE)
        }
    return(getAccuracy(fc, series))
    }
fitModel <- cmpfun(fitModel, options = list(optimize = 3))

fitThiefs <- function(series, lambda = FALSE){
    models <- c("a", "e", "f", "n", "s", "t", "z")
    results <- data.frame(a = NA, e = NA, f = NA, n = NA, s = NA, t = NA, z = NA)
    h <- length(series$xx)
    if(frequency(series$x) == 1){
        return(results)
        }
    aMod <- thief(series$x, h = h, usemodel = "arima")
    results$a[1] <- getAccuracy(aMod, series)
    eMod <- thief(series$x, h = h, usemodel = "ets")
    results$e[1] <- getAccuracy(eMod, series)
    fMod <- thief(series$x, h = h, usemodel = "theta")
    results$f[1] <- getAccuracy(fMod, series)
    nMod <- try({thief(series$x, h = h,
                      forecastfunction = function(y, h, ...) forecast(nnetar(y), h = h))},
                silent = TRUE)
    if(class(nMod) == "forecast"){
        results$n[1] <- getAccuracy(nMod, series)
        }
    #s <- thief(series$x, h = length(x$xx), usemodel = "arima")
    #s <- as.numeric(accuracy(a, x = series$xx)["Test set", "MASE"])
    tMod <- try({thief(series$x, h = h,
                      forecastfunction = function(y, h, ...) forecast(tbats(y), h = h))},
                silent = TRUE)
    if(class(tMod) == "forecast"){
        results$t[1] <- getAccuracy(tMod, series)
        }
    zMod <- thief(series$x, h = h, usemodel = "snaive")
    results$z[1] <- getAccuracy(zMod, series)
    #return(data.frame(results))
    return(results)
    }
fitThiefs <- cmpfun(fitThiefs, options = list(optimize = 3))

prophetDaily <- function(x){
    dates <- seq(from = as.Date("2000-01-01"), by = "day", length.out = length(x))
    df <- data.frame(ds = dates, y = x)
    m <- prophet(df)
    future <- make_future_dataframe(m, periods = 365)
    fc <- predict(m, future)
    return(fc)
    }
prophetDaily <- cmpfun(prophetDaily, options = list(optimize = 3))
