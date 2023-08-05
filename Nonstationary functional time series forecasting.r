###################################################
# Nonstationary functional time series forecasting
###################################################

###############################
## load and install R package
###############################

require(demography)
require(ftsa)
require(imputeTS)
require(ggplot2)
require(dplyr)
require(tidyr)
require(stringr)
require(vars)
require(doMC)

#############
## read data
#############

SWE_demo = extract.ages(read.demogdata(file = "SWE_mort.R", popfile = "SWE_expo.R", type = "mortality", label = "SWE"), 0:95)

# set those rates > 1 to 1

SWE_demo$rate$female = replace(SWE_demo$rate$female, which(SWE_demo$rate$female > 1), NA)
SWE_demo$rate$male   = replace(SWE_demo$rate$male,   which(SWE_demo$rate$male > 1),   NA)

####################
## smooth demogdata
####################

# SWE_demo_smooth <- smooth.demogdata(SWE_demo)

# set those 0 values to NA

SWE_demo$rate$female = replace(SWE_demo$rate$female, which(SWE_demo$rate$female == 0), NA)
SWE_demo$rate$male   = replace(SWE_demo$rate$male,   which(SWE_demo$rate$male == 0), NA)

# na.interp those NAs

n_age = nrow(SWE_demo$rate$female)  # 96
n_year = ncol(SWE_demo$rate$female) # 272

# linear interpolation for each age over years

SWE_demo_rate_female = SWE_demo_rate_male = matrix(NA, n_age, n_year)
for(ik in 1:n_age)
{
    SWE_demo_rate_female[ik,] = na_interpolation(SWE_demo$rate$female[ik,])
    SWE_demo_rate_male[ik,]   = na_interpolation(SWE_demo$rate$male[ik,])
    rm(ik)
}

SWE_demo_female_log = log(SWE_demo_rate_female, base = 10)
SWE_demo_male_log   = log(SWE_demo_rate_male, base = 10)
colnames(SWE_demo_female_log) = colnames(SWE_demo_male_log) = SWE_demo$year
rownames(SWE_demo_female_log) = rownames(SWE_demo_male_log) = SWE_demo$age

# check

range(SWE_demo_rate_female) # 1.6e-5 0.827922
range(SWE_demo_rate_male)   # 1.5e-5 0.9797953
which(!is.finite(SWE_demo_female_log)) # NA
which(!is.finite(SWE_demo_male_log))   # NA

########################################################
# Geometrically decaying weighted (long-run covariance)
# The procedure of Isreal Martinez
########################################################

# eigenvalue ratio criterion

select_K <- function(tau, eigenvalue)
{
    k_max = length(eigenvalue)
    k_all = rep(0, k_max-1)
    for(k in 1:(k_max-1))
    {
      k_all[k] = (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
    }
    K_hat = which.min(k_all)
    return(K_hat)
}


GW_LRC_nonstationary <- function(beta, dat, fh, fmethod = c("ets", "arima"))
{
    fmethod = match.arg(fmethod)
    n = ncol(dat)
    n_dim = nrow(dat)
    
    # equal weight or geometrically decaying weight
    
    if(is.null(beta))
    {
        q = rep(1/n, n)
    }
    else
    {
        q = beta * (1 - beta)^(0:(n - 1))
    }
    
    # compute mean
    
    # my = apply(dat, 1, weighted.mean, w = rev(q))
    my = rowMeans(dat)
    
    # de-center data by subtracting the mean term
    
    new_dat = sweep(x = t(dat), MARGIN = 2, STATS = my, FUN = "-")
    
    # first-order differencing of de-centered data
    
    new_dat_diff = matrix(NA, (n-1), nrow(dat))
    for(ik in 1:(n-1))
    {
        new_dat_diff[ik,] = new_dat[(ik+1),] - new_dat[ik,]
    }
    
    # equal weighted or geometrically decaying weighted de-centered data
    
    wq = diag(rev(q)[-n])
    new_dat2 = wq %*% new_dat_diff * (n-1)
    
    # compute long-run covariance for differenced data
    
    LRC_est = long_run_covariance_estimation(dat = t(new_dat2))
    eigen_decomp = eigen(LRC_est, symmetric = TRUE)
    
    # determine ncomp by the ratio method
    
    lambda_val = eigen_decomp$values[which(eigen_decomp$values > 0)]
    ncomp = select_K(tau = 1/log(max(lambda_val[1], length(lambda_val))), eigenvalue = lambda_val)
    rm(lambda_val); rm(LRC_est)
    
    # compute the basis function and their scores
    
    LRC_basis = matrix(eigen_decomp$vectors[,1:ncomp], ncol = ncomp)
    LRC_score = new_dat %*% LRC_basis
    
    LRC_recon = LRC_basis %*% t(LRC_score)
    LRC_resi = t(new_dat) - LRC_recon
    
    LRC_score_fore = matrix(NA, fh, ncomp)
    if(fmethod == "ets")
    {
        for(ik in 1:ncomp)
        {
            LRC_score_fore[,ik] = forecast(ets(as.numeric(LRC_score[,ik])), h = fh)$mean
        }
    }
    else if(fmethod == "arima")
    {
        for(ik in 1:ncomp)
        {
            LRC_score_fore[,ik] = forecast(auto.arima(as.numeric(LRC_score[,ik])), h = fh)$mean
        }
    }
    else
    {
        warning("Please choose a correct forecasting method.")
    }
    
    # compute long-run covariance for the residual functions
    
    LRC_resi_est = long_run_covariance_estimation(dat = LRC_resi)
    eigen_decomp_resi = eigen(LRC_resi_est, symmetric = TRUE)
    
    # determine ncomp by the ratio method
    
    lambda_val = eigen_decomp_resi$values[which(eigen_decomp_resi$values > 0)]
    ncomp_resi = select_K(tau = 1/log(max(lambda_val[1], length(lambda_val))), eigenvalue = lambda_val)
    rm(lambda_val); rm(LRC_resi_est)
    
    # compute the additional PCs and their scores
    
    LRC_resi_basis = as.matrix(eigen_decomp_resi$vectors[,1:ncomp_resi])
    LRC_resi_score = new_dat %*% LRC_resi_basis
    
    LRC_resi_score_fore = matrix(NA, fh, ncomp_resi)
    if(fmethod == "ets")
    {
        for(ik in 1:ncomp_resi)
        {
            LRC_resi_score_fore[,ik] = forecast(ets(as.numeric(LRC_resi_score[,ik])), h = fh)$mean
        }
    }
    else if(fmethod == "arima")
    {
        for(ik in 1:ncomp_resi)
        {
            LRC_resi_score_fore[,ik] = forecast(auto.arima(as.numeric(LRC_resi_score[,ik])), h = fh)$mean
        }
    }
    else
    {
        warning("Please choose a correct forecasting method.")
    }
    
    LRC_fore_curve = LRC_basis %*% matrix(LRC_score_fore, nrow = ncomp) + LRC_resi_basis %*% matrix(LRC_resi_score_fore, nrow = ncomp_resi)
    LRC_fore_curve_add_mean = LRC_fore_curve + matrix(rep(my, fh),  length(my), fh, byrow = FALSE)
    return(LRC_fore_curve_add_mean[,fh])
}

############################################################
# Selection of the optimal parameter for point forecasts
# beta_val: beta parameter
# criterion: "LRC_nonstationary"
# data: data set
# horizon: forecast horizon
# validation_year: year for validation
# err: forecast error
# forecasting_method: Forecasting_method = "ETS" or "ARIMA"
############################################################

parameter_optim <- function(beta_val, data, horizon, validation_year, err, forecasting_method)
{
    grid = nrow(data)
    n_validation_year = length(validation_year)
    validation_fore = matrix(NA, grid, (n_validation_year + 1 - horizon))
    for(iwk in 1:(n_validation_year + 1 - horizon))
    {
        validation_fore[,iwk] = GW_LRC_nonstationary(beta = beta_val, dat = data[,1:(min(validation_year) - 2 + iwk)],
                                                     fh = horizon, fmethod = forecasting_method)
    }
    if(err == "RMSPE")
    {
        err = ftsa:::rmspe(forecast = validation_fore, true = data[,(min(validation_year) + horizon - 1):max(validation_year)])
    }
    else if(err == "MAPE")
    {
        err = ftsa:::mape(forecast = validation_fore, true = data[,(min(validation_year) + horizon - 1):max(validation_year)])
    }
    else
    {
        warning("error is not correctly specified.")
    }
    return(err)
}


### RMSPE

## female

# ETS 

optim_mat_RMSPE_SWE_female_LRC_nonstationary = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim, lower = 0, upper = 1,
                   data = SWE_demo_female_log, horizon = iwjk, 
                   validation_year = SWE_validation_year, err = "RMSPE", forecasting_method = "ets")
    optim_mat_RMSPE_SWE_female_LRC_nonstationary[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}

# ARIMA

optim_mat_RMSPE_SWE_female_LRC_nonstationary_ARIMA = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim, lower = 0, upper = 1,
                   data = SWE_demo_female_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "RMSPE", forecasting_method = "arima")
    optim_mat_RMSPE_SWE_female_LRC_nonstationary_ARIMA[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}

## male

# ETS

optim_mat_RMSPE_SWE_male_LRC_nonstationary = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim, lower = 0, upper = 1,
                   data = SWE_demo_male_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "RMSPE", forecasting_method = "ets")
    optim_mat_RMSPE_SWE_male_LRC_nonstationary[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}

# ARIMA

optim_mat_RMSPE_SWE_male_LRC_nonstationary_ARIMA = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim, lower = 0, upper = 1,
                   data = SWE_demo_male_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "RMSPE", forecasting_method = "arima")
    optim_mat_RMSPE_SWE_male_LRC_nonstationary_ARIMA[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}


### MAPE

## female

# ETS

optim_mat_MAPE_SWE_female_LRC_nonstationary = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim, lower = 0, upper = 1,
                   data = SWE_demo_female_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "MAPE", forecasting_method = "ets")
    optim_mat_MAPE_SWE_female_LRC_nonstationary[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}

# ARIMA

optim_mat_MAPE_SWE_female_LRC_nonstationary_ARIMA = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim, lower = 0, upper = 1,
                   data = SWE_demo_female_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "MAPE", forecasting_method = "arima")
    optim_mat_MAPE_SWE_female_LRC_nonstationary_ARIMA[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}  

## male

# ETS

optim_mat_MAPE_SWE_male_LRC_nonstationary = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim, lower = 0, upper = 1,
                   data = SWE_demo_male_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "MAPE", forecasting_method = "ets")
    optim_mat_MAPE_SWE_male_LRC_nonstationary[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}

# ARIMA

optim_mat_MAPE_SWE_male_LRC_nonstationary_ARIMA = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim, lower = 0, upper = 1,
                   data = SWE_demo_male_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "MAPE", forecasting_method = "arima")
    optim_mat_MAPE_SWE_male_LRC_nonstationary_ARIMA[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}



############################################################
# Selection of the optimal parameter for interval forecasts
# beta_val: beta parameter
# data: data set
# horizon: forecast horizon
# validation_year: year for validation
# err: forecast error
# no_boot: number of bootstrap samples
# PI_level: nominal coverage probability
# forecasting_method: Forecasting_method = "ETS" or "ARIMA"
############################################################

parameter_optim_interval <- function(beta_val, data, horizon, validation_year, err, no_boot, 
                                     PI_level, forecasting_method)
{
    grid = nrow(data)
    n_validation_year = length(validation_year)
    test_fore_lb = test_fore_ub = matrix(NA, grid, (n_validation_year + 1 - horizon))
    for(iwk in 1:(n_validation_year + 1 - horizon))
    {
        dum = GW_LRC_nonstationary_boot(beta = beta_val, dat = data[,1:(min(validation_year) - 2 + iwk)], 
                                        fh = horizon, B = no_boot, level = PI_level, fmethod = forecasting_method)
        test_fore_lb[,iwk] = dum$lb[,horizon]
        test_fore_ub[,iwk] = dum$ub[,horizon]
        rm(iwk); rm(dum)
    }
    rownames(test_fore_lb) = rownames(test_fore_ub) = 0:(grid - 1)
    dum = interval_score(holdout = data[,(min(validation_year) + horizon - 1):max(validation_year)], 
                         lb = test_fore_lb, ub = test_fore_ub, alpha = (100 - PI_level)/100)
    if(err == "CPD")
    {
        return(dum[3])
    }
    else if(err == "score")
    {
        return(dum[1])
    }
    else
    {
        warning("err should be interval score or CPD.")
    }
}

# ETS

optim_mat_score_SWE_female_ETS = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim_interval, lower = 0, upper = 1, 
                   data = SWE_demo_female_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "score", no_boot = 1000, PI_level = 80,
                   forecasting_method = "ets")
    optim_mat_score_SWE_female_ETS[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}
rownames(optim_mat_score_SWE_female_ETS) = 1:30
colnames(optim_mat_score_SWE_female_ETS) = c("Parameter", "Objective")

# ARIMA

optim_mat_score_SWE_female_ARIMA = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim_interval, lower = 0, upper = 1, 
                   data = SWE_demo_female_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "score", no_boot = 1000, PI_level = 80,
                   forecasting_method = "arima")
    optim_mat_score_SWE_female_ARIMA[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}
rownames(optim_mat_score_SWE_female_ARIMA) = 1:30
colnames(optim_mat_score_SWE_female_ARIMA) = c("Parameter", "Objective")

# ETS

optim_mat_score_SWE_male_ETS = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim_interval, lower = 0, upper = 1, 
                   data = SWE_demo_male_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "score", no_boot = 1000, PI_level = 80,
                   forecasting_method = "ets")
    optim_mat_score_SWE_male_ETS[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}
rownames(optim_mat_score_SWE_male_ETS) = 1:30
colnames(optim_mat_score_SWE_male_ETS) = c("Parameter", "Objective")

# ARIMA

optim_mat_score_SWE_male_ARIMA = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim_interval, lower = 0, upper = 1, 
                   data = SWE_demo_male_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "score", no_boot = 1000, PI_level = 80,
                   forecasting_method = "arima")
    optim_mat_score_SWE_male_ARIMA[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}
rownames(optim_mat_score_SWE_male_ARIMA) = 1:30
colnames(optim_mat_score_SWE_male_ARIMA) = c("Parameter", "Objective")

######
# CPD
######

## female

# ETS

optim_mat_CPD_SWE_female_ETS = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim_interval, lower = 0, upper = 1,
                   data = SWE_demo_female_log, horizon = iwjk, 
                   validation_year = SWE_validation_year, err = "CPD", no_boot = 1000, PI_level = 80, 
                   forecasting_method = "ets")
    optim_mat_CPD_SWE_female_ETS[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}
rownames(optim_mat_CPD_SWE_female_ETS) = 1:30
colnames(optim_mat_CPD_SWE_female_ETS) = c("Parameter", "Objective")

# ARIMA

optim_mat_CPD_SWE_female_ARIMA = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim_interval, lower = 0, upper = 1, 
                   data = SWE_demo_female_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "CPD", no_boot = 1000, PI_level = 80,
                   forecasting_method = "arima")
    optim_mat_CPD_SWE_female_ARIMA[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}
rownames(optim_mat_CPD_SWE_female_ARIMA) = 1:30
colnames(optim_mat_CPD_SWE_female_ARIMA) = c("Parameter", "Objective")

## male

# ETS

optim_mat_CPD_SWE_male_ETS = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim_interval, lower = 0, upper = 1,
                   data = SWE_demo_male_log, horizon = iwjk, 
                   validation_year = SWE_validation_year, err = "CPD", no_boot = 1000, PI_level = 80, 
                   forecasting_method = "ets")
    optim_mat_CPD_SWE_male_ETS[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}
rownames(optim_mat_CPD_SWE_male_ETS) = 1:30
colnames(optim_mat_CPD_SWE_male_ETS) = c("Parameter", "Objective")

# ARIMA

optim_mat_CPD_SWE_male_ARIMA = matrix(NA, 30, 2)
for(iwjk in 1:30)
{
    dum = optimize(f = parameter_optim_interval, lower = 0, upper = 1, 
                   data = SWE_demo_male_log, horizon = iwjk,
                   validation_year = SWE_validation_year, err = "CPD", no_boot = 1000, PI_level = 80,
                   forecasting_method = "arima")
    optim_mat_CPD_SWE_male_ARIMA[iwjk,] = c(dum$minimum, dum$objective)
    print(iwjk); rm(iwjk); rm(dum)
}
rownames(optim_mat_CPD_SWE_male_ARIMA) = 1:30
colnames(optim_mat_CPD_SWE_male_ARIMA) = c("Parameter", "Objective")

#####################
#####################
# Empirical Analysis
#####################
#####################

########################################
# Tuning parameters
# training_year: 80% of the full data
# validation_year: 10% of the full data
# test_year: 10% of the full data
########################################

SWE_training_year = 1:212
SWE_validation_year = 213:242
SWE_test_year = 243:272

#################
# Testing sample
#################

# beta_selected: selected tuning parameter
# criterion: dynamic, two-stage dynamic, static FPCA
# data: p by n data matrix
# horizon: forecast horizon
# test_year: data in the testing period
# err: error criterion
# forecasting_method: univariate time-series forecasting method

point_fore_eval <- function(beta_selected, data, horizon, test_year, err, forecasting_method)
{
    grid = nrow(data)
    n_test_year = length(test_year)
    test_fore = matrix(NA, grid, (n_test_year + 1 - horizon))
    for(iwk in 1:(n_test_year + 1 - horizon))
    {
        test_fore[,iwk] = GW_LRC_nonstationary(beta = beta_selected, dat = data[,1:(min(test_year) - 2 + iwk)], fh = horizon,
                                               fmethod = forecasting_method)
        print(iwk); rm(iwk)
    }
    if(err == "RMSPE")
    {
        err = ftsa:::rmspe(forecast = test_fore, true = data[,(min(test_year) + horizon - 1):max(test_year)])
    }
    else if(err == "MAPE")
    {
        err = ftsa:::mape(forecast = test_fore, true = data[,(min(test_year) + horizon - 1):max(test_year)])
    }
    else
    {
        warning("error is not correctly specified.")
    }
    return(err)
}

##################
## RMSPE (Female)
# original scale
##################

## LRC_nonstationary (beta = NULL implying equal weighting)

# ETS

RMSPE_SWE_female_LRC_nonstationary = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    RMSPE_SWE_female_LRC_nonstationary[iwjk] = point_fore_eval(beta_selected = NULL,
                                                               data = SWE_demo_female_log, horizon = iwjk,
                                                               test_year = SWE_test_year, err = "RMSPE", forecasting_method = "ets")
    print(iwjk); rm(iwjk)
}
round(mean(RMSPE_SWE_female_LRC_nonstationary), 4) # 8.6317

# ARIMA

RMSPE_SWE_female_LRC_nonstationary_ARIMA = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    RMSPE_SWE_female_LRC_nonstationary_ARIMA[iwjk] = point_fore_eval(beta_selected = NULL,
                                                               data = SWE_demo_female_log, horizon = iwjk,
                                                               test_year = SWE_test_year, err = "RMSPE", forecasting_method = "arima")
    print(iwjk); rm(iwjk)
}

round(mean(RMSPE_SWE_female_LRC_nonstationary_ARIMA), 4) # 6.719

## LRC_nonstationary (beta = selected)

# ETS

RMSPE_SWE_female_LRC_nonstationary_weighted = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    RMSPE_SWE_female_LRC_nonstationary_weighted[iwjk] = point_fore_eval(beta_selected = optim_mat_RMSPE_SWE_female_LRC_nonstationary[iwjk,1],
                                                                        data = SWE_demo_female_log,
                                                                        horizon = iwjk, test_year = SWE_test_year,
                                                                        err = "RMSPE", forecasting_method = "ets")
    print(iwjk); rm(iwjk)
}
round(mean(RMSPE_SWE_female_LRC_nonstationary_weighted), 4) # 8.1859

# ARIMA

RMSPE_SWE_female_LRC_nonstationary_weighted_ARIMA = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    RMSPE_SWE_female_LRC_nonstationary_weighted_ARIMA[iwjk] = point_fore_eval(beta_selected = optim_mat_RMSPE_SWE_female_LRC_nonstationary_ARIMA[iwjk,1],
                                                                        data = SWE_demo_female_log,
                                                                        horizon = iwjk, test_year = SWE_test_year,
                                                                        err = "RMSPE", forecasting_method = "arima")
    print(iwjk); rm(iwjk)
}

round(mean(RMSPE_SWE_female_LRC_nonstationary_weighted_ARIMA), 4) # 7.3512

###################
# summary (female)
###################

RMSPE_SWE_female_summary = data.frame(cbind(RMSPE_SWE_female_LRC_nonstationary,
                                            RMSPE_SWE_female_LRC_nonstationary_weighted,
                                            RMSPE_SWE_female_LRC_nonstationary_ARIMA,
                                            RMSPE_SWE_female_LRC_nonstationary_weighted_ARIMA))
colnames(RMSPE_SWE_female_summary) = c("LRC_nonstationary", "LRC_nonstationary_weighted",
                                       "LRC", "LRC_weighted")
rownames(RMSPE_SWE_female_summary) = 1:30


###############
## RMSPE (Male)
###############

## LRC_nonstationary (beta = NULL)

# ETS

RMSPE_SWE_male_LRC_nonstationary = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    RMSPE_SWE_male_LRC_nonstationary[iwjk] = point_fore_eval(beta_selected = NULL,
                                                             data = SWE_demo_male_log, horizon = iwjk,
                                                             test_year = SWE_test_year, err = "RMSPE", forecasting_method = "ets")
    print(iwjk); rm(iwjk)
}
round(mean(RMSPE_SWE_male_LRC_nonstationary), 4) # 10.0965

# ARIMA

RMSPE_SWE_male_LRC_nonstationary_ARIMA = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    RMSPE_SWE_male_LRC_nonstationary_ARIMA[iwjk] = point_fore_eval(beta_selected = NULL,
                                                             data = SWE_demo_male_log, horizon = iwjk,
                                                             test_year = SWE_test_year, err = "RMSPE", forecasting_method = "arima")
    print(iwjk); rm(iwjk)
}

round(mean(RMSPE_SWE_male_LRC_nonstationary_ARIMA), 4) # 8.4524

## LRC_nonstationary (beta = weighted)

# ETS

RMSPE_SWE_male_LRC_nonstationary_weighted = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    RMSPE_SWE_male_LRC_nonstationary_weighted[iwjk] = point_fore_eval(beta_selected = optim_mat_RMSPE_SWE_male_LRC_nonstationary[iwjk,1],
                                                                      data = SWE_demo_male_log,
                                                                      horizon = iwjk, test_year = SWE_test_year,
                                                                      err = "RMSPE", forecasting_method = "ets")
    print(iwjk); rm(iwjk)
}
round(mean(RMSPE_SWE_male_LRC_nonstationary_weighted), 4) # 9.0875

# ARIMA

RMSPE_SWE_male_LRC_nonstationary_weighted_ARIMA = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    RMSPE_SWE_male_LRC_nonstationary_weighted_ARIMA[iwjk] = point_fore_eval(beta_selected = optim_mat_RMSPE_SWE_male_LRC_nonstationary_ARIMA[iwjk,1],
                                                                      data = SWE_demo_male_log,
                                                                      horizon = iwjk, test_year = SWE_test_year,
                                                                      err = "RMSPE", forecasting_method = "arima")
    print(iwjk); rm(iwjk)
}

round(mean(RMSPE_SWE_male_LRC_nonstationary_weighted_ARIMA), 4) # 7.9802

# summary (female)

RMSPE_SWE_male_summary = data.frame(cbind(RMSPE_SWE_male_LRC_nonstationary, RMSPE_SWE_male_LRC_nonstationary_weighted,
                                          RMSPE_SWE_male_LRC_nonstationary_ARIMA, RMSPE_SWE_male_LRC_nonstationary_weighted_ARIMA))
colnames(RMSPE_SWE_male_summary) = c("LRC_nonstationary", "LRC_nonstationary_weighted",
                                     "LRC", "LRC_weighted")
rownames(RMSPE_SWE_male_summary) = 1:30


#################
## MAPE (Female)
#################

## LRC_nonstationary (beta = NULL)

# ETS

MAPE_SWE_female_LRC_nonstationary = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    MAPE_SWE_female_LRC_nonstationary[iwjk] = point_fore_eval(beta_selected = NULL, data = SWE_demo_female_log, horizon = iwjk,
                                                              test_year = SWE_test_year, err = "MAPE", forecasting_method = "ets")
    print(iwjk); rm(iwjk)
}
round(mean(MAPE_SWE_female_LRC_nonstationary), 4) # 6.8909

# ARIMA

MAPE_SWE_female_LRC_nonstationary_ARIMA = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    MAPE_SWE_female_LRC_nonstationary_ARIMA[iwjk] = point_fore_eval(beta_selected = NULL, data = SWE_demo_female_log, horizon = iwjk,
                                                              test_year = SWE_test_year, err = "MAPE", forecasting_method = "arima")
    print(iwjk); rm(iwjk)
}

round(mean(MAPE_SWE_female_LRC_nonstationary_ARIMA), 4) # 5.1497

## LRC_nonstationary (beta = weighted)

# ETS

MAPE_SWE_female_LRC_nonstationary_weighted = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    MAPE_SWE_female_LRC_nonstationary_weighted[iwjk] = point_fore_eval(beta_selected = optim_mat_MAPE_SWE_female_LRC_nonstationary[iwjk,1],
                                                                       data = SWE_demo_female_log,
                                                                       horizon = iwjk, test_year = SWE_test_year,
                                                                       err = "MAPE", forecasting_method = "ets")
    print(iwjk); rm(iwjk)
}
round(mean(MAPE_SWE_female_LRC_nonstationary_weighted), 4) # 6.5301

# ARIMA

MAPE_SWE_female_LRC_nonstationary_weighted_ARIMA = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    MAPE_SWE_female_LRC_nonstationary_weighted_ARIMA[iwjk] = point_fore_eval(beta_selected = optim_mat_MAPE_SWE_female_LRC_nonstationary_ARIMA[iwjk,1],
                                                                            data = SWE_demo_female_log,
                                                                            horizon = iwjk, test_year = SWE_test_year,
                                                                            err = "MAPE", forecasting_method = "arima")
    print(iwjk); rm(iwjk)
}
round(mean(MAPE_SWE_female_LRC_nonstationary_weighted_ARIMA), 4) # 5.6138

###################
# summary (female)
###################

MAPE_SWE_female_summary = data.frame(cbind(MAPE_SWE_female_LRC_nonstationary, MAPE_SWE_female_LRC_nonstationary_weighted,
                                           MAPE_SWE_female_LRC_nonstationary_ARIMA, MAPE_SWE_female_LRC_nonstationary_weighted_ARIMA))
colnames(MAPE_SWE_female_summary) = c("LRC_nonstationary", "LRC_nonstationary_weighted",
                                      "LRC", "LRC_weighted")
rownames(MAPE_SWE_female_summary) = 1:30


###############
## MAPE (Male)
###############

## LRC_nonstationary (beta = NULL)

# ETS

MAPE_SWE_male_LRC_nonstationary = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    MAPE_SWE_male_LRC_nonstationary[iwjk] = point_fore_eval(beta_selected = NULL, criterion = "LRC_nonstationary",
                                                            data = SWE_demo_male_log, horizon = iwjk,
                                                            test_year = SWE_test_year, err = "MAPE", forecasting_method = "ets")
    print(iwjk); rm(iwjk)
}
round(mean(MAPE_SWE_male_LRC_nonstationary), 4) # 8.4416

# ARIMA

MAPE_SWE_male_LRC_nonstationary_ARIMA = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    MAPE_SWE_male_LRC_nonstationary_ARIMA[iwjk] = point_fore_eval(beta_selected = NULL,
                                                            data = SWE_demo_male_log, horizon = iwjk,
                                                            test_year = SWE_test_year, err = "MAPE", forecasting_method = "arima")
    print(iwjk); rm(iwjk)
}
round(mean(MAPE_SWE_male_LRC_nonstationary_ARIMA), 4) # 6.9584

## LRC_nonstationary (beta = selected)

# ETS

MAPE_SWE_male_LRC_nonstationary_weighted = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    MAPE_SWE_male_LRC_nonstationary_weighted[iwjk] = point_fore_eval(beta_selected = optim_mat_MAPE_SWE_male_LRC_nonstationary[iwjk,1],
                                                                     data = SWE_demo_male_log,
                                                                     horizon = iwjk, test_year = SWE_test_year,
                                                                     err = "MAPE", forecasting_method = "ets")
    print(iwjk); rm(iwjk)
}
round(mean(MAPE_SWE_male_LRC_nonstationary_weighted), 4) # 7.6075

# ARIMA

MAPE_SWE_male_LRC_nonstationary_weighted_ARIMA = vector("numeric", length(SWE_test_year))
for(iwjk in 1:length(SWE_test_year))
{
    MAPE_SWE_male_LRC_nonstationary_weighted_ARIMA[iwjk] = point_fore_eval(beta_selected = optim_mat_MAPE_SWE_male_LRC_nonstationary_ARIMA[iwjk,1],
                                                                     data = SWE_demo_male_log,
                                                                     horizon = iwjk, test_year = SWE_test_year,
                                                                     err = "MAPE", forecasting_method = "arima")
    print(iwjk); rm(iwjk)
}
round(mean(MAPE_SWE_male_LRC_nonstationary_weighted_ARIMA), 4) # 6.4053

# summary (male)

MAPE_SWE_male_summary = data.frame(cbind(MAPE_SWE_male_LRC_nonstationary, MAPE_SWE_male_LRC_nonstationary_weighted,
                                         MAPE_SWE_male_LRC_nonstationary_ARIMA, MAPE_SWE_male_LRC_nonstationary_weighted_ARIMA))
colnames(MAPE_SWE_male_summary) = c("LRC_nonstationary", "LRC_nonstationary_weighted",
                                    "LRC", "LRC_weighted")
rownames(MAPE_SWE_male_summary) = 1:30

########
# Plots
########

## univariate time series plots

# females

savepdf("Fig_1a", width = 12, height = 10, toplines = 0.8)
plot(ts(log(SWE_demo$rate$female, base = 10)[1,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
     xlab = "Year", ylab = expression(paste(paste("Mortality rate (", log[10]), sep="", " scale)", sep=" ")),
     ylim = c(-4.5, 0), col = 1, lty = 1, main = "Swedish females (1751 - 2022)")
lines(ts(log(SWE_demo$rate$female, base = 10)[21,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 2, lty = 2)
lines(ts(log(SWE_demo$rate$female, base = 10)[41,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 3, lty = 3)
lines(ts(log(SWE_demo$rate$female, base = 10)[61,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 4, lty = 4)
lines(ts(log(SWE_demo$rate$female, base = 10)[81,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 5, lty = 5)
lines(ts(log(SWE_demo$rate$female, base = 10)[96,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 6, lty = 6)
legend("bottomleft", c("0", "20", "40", "60", "80", "95+"), col = 1:6, lty = 1:6, ncol = 2, cex = 0.8)
dev.off()

# males

savepdf("Fig_1b", width = 12, height = 10, toplines = 0.8)
plot(ts(log(SWE_demo$rate$male, base = 10)[1,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
     xlab = "Year", ylab = expression(paste(paste("Mortality rate (", log[10]), sep="", " scale)", sep=" ")),
     ylim = c(-4.5, 0), col = 1, lty = 1, main = "Swedish males (1751 - 2022)")
lines(ts(log(SWE_demo$rate$male, base = 10)[21,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 2, lty = 2)
lines(ts(log(SWE_demo$rate$male, base = 10)[41,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 3, lty = 3)
lines(ts(log(SWE_demo$rate$male, base = 10)[61,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 4, lty = 4)
lines(ts(log(SWE_demo$rate$male, base = 10)[81,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 5, lty = 5)
lines(ts(log(SWE_demo$rate$male, base = 10)[96,], start = min(SWE_demo$year), end = max(SWE_demo$year)),
      col = 6, lty = 6)
dev.off()

# functional time series plots

savepdf("Fig_1c", width = 12, height = 10, toplines = 0.8)
ages = 0:95
plot(fts(ages, SWE_demo_female_log), xlab = "Age",
         ylab = expression(paste(paste("Mortality rate (", log[10]), sep="", " scale)", sep=" ")),
     main = "")
dev.off()

savepdf("Fig_1d", width = 12, height = 10, toplines = 0.8)
plot(fts(ages, SWE_demo_male_log), xlab = "Age",
     ylab = "",
     main = "")
dev.off()


## Plots of forecasts

savepdf("Fig_3a", width = 12, height = 10)
RMSPE_SWE_female_summary %>%
  mutate(x=1:30) %>%
  pivot_longer(c(LRC_nonstationary, LRC_nonstationary_weighted, LRC, LRC_weighted)) %>%
  mutate(color = case_when(name %in% c("LRC_nonstationary", "LRC_nonstationary_weighted") ~ "ETS",
                            name %in% c("LRC", "LRC_weighted") ~ "ARIMA"),
         lty = case_when(str_ends(name, "_weighted") ~ "Weighted", TRUE ~ "Standard")) %>%
  ggplot(aes(x, y = value, color = color, lty = lty)) +
  geom_line() +
  xlab('') +
  ylab('RMSPE') +
  ylim(4.5, 16.5) +
  labs(title='Female') +
  theme(legend.position = c(0.12,0.72))
dev.off()


savepdf("Fig_3b", width = 12, height = 10)
RMSPE_SWE_male_summary %>%
  mutate(x=1:30) %>%
  pivot_longer(c(LRC_nonstationary, LRC_nonstationary_weighted, LRC, LRC_weighted)) %>%
  mutate(color = case_when(name %in% c("LRC_nonstationary", "LRC_nonstationary_weighted") ~ "ETS",
                           name %in% c("LRC", "LRC_weighted") ~ "ARIMA"),
         lty = case_when(str_ends(name, "_weighted") ~ "Weighted", TRUE ~ "Standard")) %>%
  ggplot(aes(x, y = value, color = color, lty = lty)) +
  geom_line() +
  xlab('') +
  ylab('') +
  ylim(4.5, 16.5) +
  labs(title='Male')+
  theme(legend.position = "none")
dev.off()


savepdf("Fig_3c", width = 12, height = 10)
MAPE_SWE_female_summary %>%
  mutate(x=1:30) %>%
  pivot_longer(c(LRC_nonstationary, LRC_nonstationary_weighted, LRC, LRC_weighted)) %>%
  mutate(color = case_when(name %in% c("LRC_nonstationary", "LRC_nonstationary_weighted") ~ "ETS",
                           name %in% c("LRC", "LRC_weighted") ~ "ARIMA"),
         lty = case_when(str_ends(name, "_weighted") ~ "Weighted", TRUE ~ "Standard")) %>%
  ggplot(aes(x, y = value, color = color, lty = lty)) +
  geom_line() +
  xlab('Forecast horizon') +
  ylab('MAPE') +
  ylim(2.5, 14) +
  theme(legend.position = "none")
dev.off()


savepdf("Fig_3d", width = 12, height = 10)
MAPE_SWE_male_summary %>%
  mutate(x=1:30) %>%
  pivot_longer(c(LRC_nonstationary, LRC_nonstationary_weighted, LRC, LRC_weighted)) %>%
  mutate(color = case_when(name %in% c("LRC_nonstationary", "LRC_nonstationary_weighted") ~ "ETS",
                           name %in% c("LRC", "LRC_weighted") ~ "ARIMA"),
         lty = case_when(str_ends(name, "_weighted") ~ "Weighted", TRUE ~ "Standard")) %>%
  ggplot(aes(x, y = value, color = color, lty = lty)) +
  geom_line() +
  xlab('Forecast horizon') +
  ylab('')+
  ylim(2.5, 14) +
  theme(legend.position = "none")
dev.off()







