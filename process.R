library(data.table)
library(forecastHybrid)

# Load all the data files
files <- dir(pattern = "*.csv")

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
    type <- determineType()
    dat$V1 <- NULL
    horizon <- getHorizon()
    dat = dat[, lapply(.SD, as.numeric)]
    return(dat)
    }
