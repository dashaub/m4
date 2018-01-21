library(forecastHybrid)
library(compiler)

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

fitThiefs <- function(series, lambda = FALSE){
    models <- c("a", "e", "f", "n", "s", "t", "z")
    h <- length(series$xx)
    if(frequency(series$x) == 1){
        nan <- rep(NA, 7)
        names(nan) <- models
        return(nan)
        }
       results <- list()
       a <- thief(series$x, h = h, usemodel = "arima")
       a <- as.numeric(accuracy(a, x = series$xx)["Test set", "MASE"])
       e <- thief(series$x, h = h, usemodel = "ets")
       e <- as.numeric(accuracy(e, x = series$xx)["Test set", "MASE"])
       f <- thief(series$x, h = h, usemodel = "theta")
       f <- as.numeric(accuracy(f, x = series$xx)["Test set", "MASE"])
       n <- thief(series$x, h = h,
                  forecastfunction = function(y, h, ...) forecast(nnetar(y), h = h))
       n <- as.numeric(accuracy(n, x = series$xx)["Test set", "MASE"])
       #s <- thief(series$x, h = length(x$xx), usemodel = "arima")
       #s <- as.numeric(accuracy(a, x = series$xx)["Test set", "MASE"])
       s <- NA
       t <- thief(series$x, h = h,
                  forecastfunction = function(y, h, ...) forecast(tbats(y), h = h))
       t <- as.numeric(accuracy(t, x = series$xx)["Test set", "MASE"])
       z <- thief(series$x, h = h, usemodel = "snaive")
       z <- as.numeric(accuracy(z, x = series$xx)["Test set", "MASE"])
       results <- c(a, e, f, n, s, t, z)
       names(results) <- models
    return(results)
    }
fitThiefs <- cmpfun(fitThiefs, options = list(optimize = 3))
