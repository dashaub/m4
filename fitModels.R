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
