


determineType <- function(df){
    # For each data file, determine if it is daily, hourly, weekly, etc
    return(substr(head(df$V1, 1), start = 0, 1))
    }

getHorizon <- function(x){
    # Determine the forecast horizon
    horizon <- c("HOURLY" = 48, "DAILY" = 14, "WEEKLY" = 13, "MONTHLY" = 18, "YEARLY" = 6)
    return(horizon[x])
    }
getHorizonFromFrequency <- function(x){
    freq <- as.character(frequency(x))
    switch(freq, "24" = 48, "365" = 14, "52" = 13, "12" = 18, "1" = 6)
    }

prepareM <- function(data, period, type, names){
    seriesList <- list()
    horizon <- getHorizon(period)
    count <- 1
    for(series in data){
        endSlice = length(series) - horizon
        trainSeries <- subset(series, end = endSlice)
        testSeries <- subset(series, start = endSlice + 1)
        seriesList[[seriesNames[count]]]$x <- trainSeries
        seriesList[[seriesNames[count]]]$xx <- testSeries
        seriesList[[seriesNames[count]]]$period <- period
        seriesList[[seriesNames[count]]]$type <- "MICRO"
        count <- count + 1
        }
    return(seriesList)
    }

prepareDatasets <- function(){
    # Prepare data
    datasets <- list(nn3, nn5, nngc1, gefcom2012_load, gefcom2012_temp, gefcom2012_wp)
    seriesNames <- names(nn3)
    NN3 <- prepareM(data = nn3, period = "MONTHLY", type = "OTHER", names = seriesNames)
    seriesNames <- names(nn5)
    NN5 <- prepareM(data = nn5, period = "DAILY", type = "OTHER", names = seriesNames)
    # nngc1 placeholder
    

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
