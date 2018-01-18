featuresHelper <- function(x){
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

    #ACF features
    pacf_fit <- pacf(x, plot = FALSE, lag.max = 7)
    acf_fit <- acf(x, plot = FALSE, lag.max = 7)
    acf_sig <- sum(abs(acf_fit$acf) > 2 / sqrt(len))
    pacf_sig <- sum(abs(pacf_fit$acf) > 2 / sqrt(len))
    acf_sum <- sum(abs(acf_fit$acf))
    pacf_sum <- sum(abs(pacf_fit$acf))

    # BDS
    bds_test <- bds.test(x)

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
                    ar_7 = ar_fit[7],
                    pacf_1 = pacf_fit$acf[1],
                    pacf_2 = pacf_fit$acf[2],
                    pacf_3 = pacf_fit$acf[3],
                    pacf_4 = pacf_fit$acf[4],
                    pacf_5 = pacf_fit$acf[5],
                    pacf_6 = pacf_fit$acf[6],
                    pacf_7 = pacf_fit$acf[7],
                    # Ignore cov(x_t, x_t)
                    acf_1 = acf_fit$acf[2],
                    acf_2 = acf_fit$acf[3],
                    acf_3 = acf_fit$acf[4],
                    acf_4 = acf_fit$acf[5],
                    acf_5 = acf_fit$acf[6],
                    acf_6 = acf_fit$acf[7],
                    acf_7 = acf_fit$acf[8],
                    acf_sig = acf_sig,
                    pacf_sig = pacf_sig,
                    acf_sum = acf_sum,
                    pacf_sum = pacf_sum,
                    sharpe = sharpe(scaled),
                    sterling = sterling(scaled),
                    runs_pos <- runs.test(factor(diffSeries >= 0))$p.value,
                    runs_pos_less <- runs.test(factor(diffSeries >= 0), alternative = "less")$p.value,
                    runs_pos_greater <- runs.test(factor(diffSeries >= 0), alternative = "greater")$p.value,
                    runs_neg <- runs.test(factor(diffSeries <= 0))$p.value,
                    runs_neg_less <- runs.test(factor(diffSeries <= 0), alternative = "less")$p.value,
                    runs_neg_greater <- runs.test(factor(diffSeries <= 0), alternative = "greater")$p.value,
                    bds_max <- max(bds_test$pvalue),
                    bds_min <- max(bds_test$pvalue),
                    bds_median <- mediain(bds_test$pvalue),
                    bds_mean <- mean(bds_test$pvalue)
                    )
    rownames(df) <- NULL
    options(warn=0)
    return(df)
    }

tsFeatures <- function(x){
    noDiff <- featuresHelper(x)
    diffSeries <- diff(x)
    diffFeatures <- featuresHelper(diffSeries)
    names(diffFeatures) <- paste0("diff_", names(diffFeatures))
    return(cbind(noDiff, diffFeatures))
    }
