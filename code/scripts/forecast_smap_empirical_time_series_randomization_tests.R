# Randomization test for correlations between rankings and species
# forecast errors in empirical time series 

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(plotly)) {install.packages("plotly"); library(plotly)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(corrplot)) {install.packages("corrplot"); library(corrplot)}
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require(ppcor)) {install.packages("ppcor"); library(ppcor)}
if(!require(rEDM)) {install.packages("rEDM"); library(rEDM)}
if(!require(reticulate)) {install.packages("reticulate"); library(reticulate)}
if(!require(zoo)) {install.packages("zoo"); library(zoo)}
if(!require(expm)) {install.packages("expm"); library(expm)}

# loading time series, s-map results, and forecast results ------------------------------ 
# proportion of data as training set
training_prop <- 0.7
# how many steps ahead forecast
steps_ahead <- 3
# time steps of perturbed points
k <- steps_ahead
# empirical data set to use
data <- "beninca_2009"
# percentile of leading eigenvalue to use
lambda1_perc <- 0
# load empirical time series
ts <- read.csv(paste("data/empirical_time_series/", data, ".csv", sep = ""), header = TRUE)
sp_names <- names(ts)[-1]
n_sp <- ncol(ts) - 1
# load s-map results
load(paste("results/forecast_analyses/sequential_smap_", data, "_",
           training_prop, "_training_", steps_ahead, "_steps_ahead",
           ".RData", sep = ""))
load(paste("results/forecast_analyses/sequential_smap_jacobians_", data, "_",
           training_prop, "_training_", steps_ahead, "_steps_ahead",
           ".RData", sep = ""))
# load forecast results
load(paste("results/forecast_analyses/sequential_forecast_", data, "_",
           training_prop, "_training_", steps_ahead, "_steps_ahead",
           ".RData", sep = ""))

# computing species expected sensitivities ------------------------------ 
start_last_train <- floor(nrow(ts) * training_prop)
# loop over all s-map jacobians
expected_sensitivities <- list()
for (l in 1:length(smap_jacobians)) {
  # initial covariance matrix of perturbations
  Sigma_initial <- diag(rep(1, n_sp))
  # exponential of k*J
  M <- expm(k * smap_jacobians[[l]])
  # final covariance matrix of perturbations
  Sigma_final <- M %*% Sigma_initial %*% t(M) 
  # expected sensitivities
  expected_sensitivities[[l]] <- diag(Sigma_final)
  # rescaling expected sensitivities to sum 1
  expected_sensitivities[[l]] <- expected_sensitivities[[l]] / sum(expected_sensitivities[[l]])
}
# create data frame with expected sensitivity results
exp_sens_df <- data.frame(matrix(unlist(expected_sensitivities), nrow = length(expected_sensitivities), 
                                 byrow = TRUE))
names(exp_sens_df) <- names(ts)[-1]
# last training point to start with
start_last_train <- which(ts$time == smap_df$time[1])
# first test set point
start_first_test <- which(ts$time == smap_df$time[2])

# computing forecast errors sequentially ------------------------------ 
# total number of forecasts
n_forecasts <- nrow(smap_df)
# objects to store results
rmse_ratio <- list()
# loop to compute RMSE of sequential forecasts
for (j in 1:n_forecasts) {
  # forecasts for time smap$time[j]
  forecast <- forecast_df[forecast_df$last_train_time == smap_df$time[j], -c(1, 2)]
  # just test set
  ts_test <- ts[(start_first_test + j - 1):(start_first_test + j + steps_ahead - 2), -1]
  # test set standard deviation
  test_sd <- apply(ts_test, 2, sd)
  # add a small value when sd is zero
  test_sd[test_sd == 0] <- 0.001
  # forecast RMSE
  forecast_rmse <- apply((ts_test - forecast)^2, 2, function(x) sqrt(mean(x)))
  # naive forecasts
  naive_forecast <- matrix(rep(as.numeric(ts[ts$time == smap_df$time[j], -1]), each = steps_ahead), 
                           nrow = steps_ahead,
                           ncol = n_sp, byrow = FALSE)
  # add a small value when naive forecast is zero
  naive_forecast[naive_forecast == 0] <- 0.001
  # naive forecast RMSE
  naive_rmse <- apply((ts_test - naive_forecast)^2, 2, function(x) sqrt(mean(x)))
  # standardizing forecast RMSE
  rmse_ratio[[j]] <- forecast_rmse / naive_rmse
}
# transforming list to data frame
rmse_ratio_df <- as.data.frame(do.call("rbind", rmse_ratio))

# compute correlation between species forecast errors and different rankings ------------------------------ 
# obtain subset of data frames
curr_ts <- ts[start_last_train:(nrow(ts) - steps_ahead), ]
curr_ts <- curr_ts[smap_df$lambda1 > quantile(smap_df$lambda1, probs = lambda1_perc), ]
curr_rmse_ratio_df <- rmse_ratio_df[smap_df$lambda1 > quantile(smap_df$lambda1, probs = lambda1_perc), ]
curr_smap_df <- smap_df[smap_df$lambda1 > quantile(smap_df$lambda1, probs = lambda1_perc), ]
curr_exp_sens_df <- exp_sens_df[smap_df$lambda1 > quantile(smap_df$lambda1, probs = lambda1_perc), ]
# vectors to store correlations
cor_rmse_exp_sens <- c()
cor_rmse_eigenvector <- c()
cor_rmse_rate_change <- c()
cor_rmse_abundance <- c()
# compute correlations
for (m in 1:nrow(curr_rmse_ratio_df)) {
  # correlation between RMSE and species expected sensitivities
  cor_rmse_exp_sens[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                              as.numeric(curr_exp_sens_df[m, ]),
                              method = "spearman")
  # correlation between RMSE and species eigenvector alignment
  cor_rmse_eigenvector[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                                 abs(as.numeric(curr_smap_df[m, 2:(n_sp+1)])),
                                 method = "spearman")
  # correlation between RMSE and species rate of change
  rate_change <- as.numeric(abs((curr_ts[m + 1, 2:(n_sp+1)] - curr_ts[m, 2:(n_sp+1)]) / 
                                  curr_ts[m, 2:(n_sp+1)]))
  rate_change[is.na(rate_change)] <- 0
  cor_rmse_rate_change[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), rate_change,
                                 method = "spearman")
  # correlation between RMSE and species abundance
  cor_rmse_abundance[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                               -as.numeric(curr_ts[m, 2:(n_sp+1)]),
                               method = "spearman")
}
# mean correlation values
mean_cor_rmse_exp_sens <- mean(cor_rmse_exp_sens, na.rm = TRUE)
mean_cor_rmse_eigenvector <- mean(cor_rmse_eigenvector, na.rm = TRUE)
mean_cor_rmse_rate_change <- mean(cor_rmse_rate_change, na.rm = TRUE)
mean_cor_rmse_abundance <- mean(cor_rmse_abundance, na.rm = TRUE)

# compute correlation between species forecast errors and randomized rankings ------------------------------ 
mean_random_cor_rmse_exp_sens <- c()
mean_random_cor_rmse_eigenvector <- c()
mean_random_cor_rmse_rate_change <- c()
mean_random_cor_rmse_abundance <- c()
for (i in 1:1000) {
  print(i)
  # vectors to store correlations
  cor_rmse_exp_sens <- c()
  cor_rmse_eigenvector <- c()
  cor_rmse_rate_change <- c()
  cor_rmse_abundance <- c()
  # compute correlations
  for (m in 1:nrow(curr_rmse_ratio_df)) {
    # correlation between RMSE and species expected sensitivities
    cor_rmse_exp_sens[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                                sample(as.numeric(curr_exp_sens_df[m, ])),
                                method = "spearman")
    # correlation between RMSE and species eigenvector alignment
    cor_rmse_eigenvector[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                                   sample(abs(as.numeric(curr_smap_df[m, 2:(n_sp+1)]))),
                                   method = "spearman")
    # correlation between RMSE and species rate of change
    rate_change <- as.numeric(abs((curr_ts[m + 1, 2:(n_sp+1)] - curr_ts[m, 2:(n_sp+1)]) / 
                                    curr_ts[m, 2:(n_sp+1)]))
    rate_change[is.na(rate_change)] <- 0
    cor_rmse_rate_change[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                                   sample(rate_change),
                                   method = "spearman")
    # correlation between RMSE and species abundance
    cor_rmse_abundance[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                                 sample(-as.numeric(curr_ts[m, 2:(n_sp+1)])),
                                 method = "spearman")
  }
  mean_random_cor_rmse_exp_sens[i] <- mean(cor_rmse_exp_sens, na.rm = TRUE)
  mean_random_cor_rmse_eigenvector[i] <- mean(cor_rmse_eigenvector, na.rm = TRUE)
  mean_random_cor_rmse_rate_change[i] <- mean(cor_rmse_rate_change, na.rm = TRUE)
  mean_random_cor_rmse_abundance[i] <- mean(cor_rmse_abundance, na.rm = TRUE)
}
# p values
mean(mean_random_cor_rmse_exp_sens > mean_cor_rmse_exp_sens)
mean(mean_random_cor_rmse_eigenvector > mean_cor_rmse_eigenvector)
mean(mean_random_cor_rmse_rate_change > mean_cor_rmse_rate_change)
mean(mean_random_cor_rmse_abundance > mean_cor_rmse_abundance)
