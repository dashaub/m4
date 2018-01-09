


determineType <- function(df){
    # For each data file, determine if it is daily, hourly, weekly, etc
    return(substr(head(df$V1, 1), start = 0, 1))
    }

getHorizon <- function(x){
    # Determine the forecast horizon
    horizon <- c("HOURLY" = 48, "DAILY" = 14, "WEEKLY" = 13, "MONTHLY" = 18, "YEARLY" = 6)
    return(horizon[x])
    }

prepareDatasets(){
    # Prepare data
    datasets <- list(nn3, nn5, nngc1, gefcom2012_load, gefcom2012_temp, gefcom2012_wp)
    NN3 <- list()
    period <- "MONTHLY"
    horizon <- getHorizon(period)
    count <- 1
    seriesNames <- names(nn3)
    for(series in nn3){
        endSlice = length(series) - horizon
        trainSeries <- subset(series, end = endSlice)
        testSeries <- subset(series, start = endSlice + 1)
        NN3[[seriesNames[count]]]$x <- trainSeries
        NN3[[seriesNames[count]]]$xx <- testSeries
        NN3[[seriesNames[count]]]$period <- period
        NN3[[seriesNames[count]]]$type <- "MICRO"
        count <- count + 1
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
