library(data.table)
library(forecastHybrid)

determineType(df){
    return substr(head(V1$, 1), 0, 1)
    }


processFile <- function(x){
    dat <- fread(x, header = TRUE)
    type <- determineType()
    dat$V1 <- NULL
    horizon <- getHorizon()
    dat = dat[, lapply(.SD, as.numeric), by=ID]
    return dat
    }
