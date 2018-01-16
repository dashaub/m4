library(caret)
library(ranger)
library(Boruta)
library(data.table)
library(forecastHybrid)
library(tseries)
library(e1071)
library(uroot)
library(entropy)
library(nortest)
library(pbapply)

library(Mcomp)
library(doMC)
library(Tcomp)


# Load all the data files
files <- dir(pattern = "*.csv")

# All possible models
#~ models <- c("a", "e", "f", "n", "s", "t")
#~ expandedGrid <- expand.grid(rep(list(models), times = length(models)))
#~ noDupes <- unique(apply(expandedGrid, 1, FUN = function(x) sort(unique(x))))
#~ allModels <- noDupes[sapply(noDupes, FUN = function(x) length(x) >= 1)]
#~ allModels <- sapply(allModels, FUN = function(x) paste0(x, collapse=""))
#~ allModels <- c(allModels, "z")
allModels <- c("a", "e", "f", "n", "s", "t", "z")

tsFeatures <- function(x){
    options(warn=-1)
    adf <- adf.test(x)
    scaled <- scale(x)
    min_x = min(scaled)
    max_x = max(scaled)
    tsFreq <- frequency(x) > 1
    minx <- min(x)

    # Log features
    logx <- log2(x)
    if(minx <= 0){
        logx <- log2(x + abs(minx) + 10^-16)
        }
    log_var <- var(logx)
    log_range <- max(logx) - min(logx)
    log_skew <- abs(skewness(logx))
    log_kurtosis <- kurtosis(logx)

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
    len = length(x)
    unique_len = length(unique(x))
    numBins = 11
    discrete = discretize(AirPassengers, numBins = numBins)
    entrpy = entropy(x, method = "ML")
    tabled <- table(x)

    # Nortest
    ad <- ad.test(scaled)
    cvm <- cvm.test(scaled)
    lillie <- lillie.test(scaled)
    pearson <- pearson.test(scaled)
    sf <- sf.test(scaled)

    # MA features
    ma3 <- ma(x, order = 3, centre = FALSE)
    ma3_mask <- !is.na(ma3)
    ma6 <- ma(x, order = 6, centre = FALSE)
    ma6_mask <- !is.na(ma6)
    ma3_dat <- data.frame(ma3 = ma3[ma3_mask], y = x[ma3_mask])
    ma6_dat <- data.frame(ma6 = ma6[ma6_mask], y = x[ma6_mask])
    ma3_reg <- lm(y ~ ., data = ma3_dat)
    ma6_reg <- lm(y ~ ., data = ma6_dat)
    summary_ma3 <- summary(ma3_reg)
    summary_ma6 <- summary(ma6_reg)

    # AR features
    a <- ar(x)
    ar_fit <- ar(x, aic = FALSE, order.max = 7)$ar

    # Differences
    diffSeries <- diff(x)
    posDiff <- length(which(diffSeries > 0)) / length(diffSeries)
    eqDiff <- length(which(diffSeries == 0)) / length(diffSeries)
    updn <- c(0, diff(sign(diffSeries)))
    ix <- which(updn != 0)
    signChange <- length(ix) / length(x)
    diff_sd <- sd(diffSeries) / abs(mean(diffSeries))
    autocor <- cor(x[-length(x)],x[-1])

    # Build lots of features
    df = data.frame(len = length(x),
                    unique_len = len,
                    unique_ratio = unique_len / len,
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
                    num_periods = length(x) / frequency(x),
                    entropy = entrpy,
                    abs_entropy = abs(entrpy),
                    # This produces NA
                    #entropy_diff = entropy(diff(x)),
                    # Boruta rejected
                    #discrete_entropy = entropy(discrete),
                    #discrete_median = sum(head(discrete, numBins / 2)) / sum(discrete),
                    entropy_mm = entropy(tabled, method = "MM"),
                    entropy_jeffreys = entropy(tabled, method = "Jeffreys"),
                    entropy_laplace = entropy(tabled, method = "Laplace"),
                    entropy_sg = entropy(tabled, method = "SG"),
                    entropy_minimax = entropy(tabled, method = "minimax"),
                    entropy_cs = entropy(tabled, method = "CS"),
                    spectral = findfrequency(x),
                    ad_statistic <- ad$statistic,
                    ad_p <- ad$p.value,
                    cvm_statistic <- cvm$statistic,
                    cvm_p <- cvm$p.value,
                    lillie_statistic <- lillie$statistic,
                    lillie_p <- lillie$p.value,
                    pearson_statistic <- pearson$statistic,
                    pearson_p <- pearson$p.value,
                    sf_statistic <- sf$statistic,
                    sf_p <- sf$p.value,
                    log_var = log_var,
                    log_range = log_range,
                    log_skew = log_skew,
                    log_kurtosis = log_kurtosis,
                    ma3_r_sq = summary_ma3$adj.r.squared,
                    ma3_r_coef = coef(summary_ma3)[2,1],
                    ma3_coef_t = abs(coef(summary_ma3)[2,3]),
                    ma6_r_sq = summary_ma6$adj.r.squared,
                    ma6_r_coef = coef(summary_ma6)[2,1],
                    ma6_coef_t = abs(coef(summary_ma6)[2,3]),
                    ar_order = a$order,
                    ar_resids = sum(abs(residuals(a)), na.rm = TRUE) / sum(abs(x)),
                    pos_diff = posDiff,
                    eq_diff = eqDiff,
                    sign_change = signChange,
                    diff_sd = diff_sd,
                    autocor <- autocor,
                    ar_1 = ar_fit[1],
                    ar_2 = ar_fit[2],
                    ar_3 = ar_fit[3],
                    ar_4 = ar_fit[4],
                    ar_5 = ar_fit[5],
                    ar_6 = ar_fit[6],
                    ar_7 = ar_fit[7]
                    )
    rownames(df) <- NULL
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
        } else if(mObj[[i]]$period == "OTHER"){
            #mObj[[i]]$period <- "DAILY"
            }
        }
    return(mObj)
    }

data(M3)
data(M1)
allData <- c(M1, M3, tourism)
if(require(tscompdata)){
	tscompdataPrepare(){
		timeseries <- list(nn3, nn5, nngc1, gefcom2012_load, gefcom2012_temp, gefcom2012_wp)
		}
	}

createMObject <- function(x, type){
	xClean <- na.interp(x)
	tsLength <- length(x)
	horizon <- getHorizonFromFrequency(xClean)
	testSet <- subset(xClean, start = tsLength - horizon + 1)
	trainSet <- subset(xClean, end = tsLength - horizon)
	returnList <- list(x = trainSet, x = testSet, h = horizon, n = length(trainSet), type = type)
	return(returnList)
	}

res <- lapply(nn3, FUN = function(x) createMObject(x, type = "MICRO"))
cleaned <- cleanM(allData)
# Shuffle data so it distributes evenly for parallel feature extraction
set.seed(31415926)
cleaned <- rclean <- sample(cleaned, size = length(cleaned), replace = FALSE)
# Remove short series
shortSeries <- sapply(cleaned, FUN = function(x) length(x$x) <= 7)
cleaned <- cleaned[!shortSeries]

# Prepare data for training

#~ cleaned <- sample(allData, 50)
#~ table(sapply(cleaned, FUN = function(x) x$type))
#~ table(sapply(cleaned, FUN = function(x) x$period))
#~ summary(sapply(cleaned, FUN = function(x) length(x$x)))

# Fit multiple models and return the MASE
fitModels <- function(x, models){
    results <- sapply(models, FUN = function(model) fitModel(x, method=model))
    df <- data.frame(matrix(results, nrow = 1))
    names(df) <- names(results)
    return(df)
    }

# Fit a single model and return the MASE
fitModel <- function(series, method){
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


# Build features and labels
cl <- makeForkCluster(8)
#features <- rbindlist(parLapplyLB(cl = cl, X = cleaned, fun = function(x) tsFeatures(x$x)))
features <- rbindlist(pblapply(X = cleaned, FUN = function(x) tsFeatures(x$x), cl = cl))
save(features, file = "features.RData", compress = "xz", compression_level = 9)

#res <- parLapplyLB(cl = cl, X = cleaned, fun = function(x) fitModels(x = x, models = allModels))
mase <- rbindlist(pblapply(X = cleaned,
                           FUN = function(x) fitModels(x = x, models = allModels),
                           cl = cl))
save(mase, file = "mase.RData")


#library(doParallel)
#registerDoParallel(8)
#registerDoMC(8)


means <- colMeans(mase, na.rm = TRUE)
sorted <- apply(mase, 1, FUN = sort)

# Create labels
labels_first <- factor(sapply(sorted, FUN = function(x) names(x[1])))
labels_second <- factor(sapply(sorted, FUN = function(x) names(x[2])))
labels_third <- factor(sapply(sorted, FUN = function(x) names(x[3])))
labels_second_worst <- factor(sapply(sorted, FUN = function(x) names(tail(x, 2)[1])))
labels_worst <- factor(sapply(sorted, FUN = function(x) names(tail(x, 1))))
save(labels_first, file = "labels_first.RData")
save(labels_second, file = "labels_second.RData")
save(labels_third, file = "labels_third.RData")
save(labels_second_worst, file = "labels_second_worst.RData")
save(labels_worst, file = "labels_worst.RData")


# Feature selection
dat <- features
dat$type <- as.character(sapply(cleaned, FUN = function(x) x$type))
set.seed(34)
b <- Boruta(x = dat, y = labels_worst, doTrace = 1)

# Train models
#registerDoMC(1)
tc <- trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "random")
numMod <- 150
seed <- 50
set.seed(seed)
rangerModFirst <- train(x = dat, y = labels_first,
                        method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModFirst, file = "rangerModFirst.RData",
     compress = "xz", compression_level = 9)
set.seed(seed)
rangerModSecond <- train(x = dat, y = labels_second,
                         method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModSecond, file = "rangerModSecond.RData",
     compress = "xz", compression_level = 9)
set.seed(seed)
rangerModThird <- train(x = dat, y = labels_third,
                        method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModThird, file = "rangerModThird.RData",
     compress = "xz", compression_level = 9)
set.seed(seed)
rangerModSecondWorst <- train(x = dat, y = labels_second_worst,
                              method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModSecondWorst, file = "rangerModSecondWorst.RData",
     compress = "xz", compression_level = 9)
set.seed(seed)
rangerModWorst <- train(x = dat, y = labels_worst,
                        method = "ranger", trControl = tc,tuneLength = numMod)
save(rangerModWorst, file = "rangerModWorst.RData",
     compress = "xz", compression_level = 9)
