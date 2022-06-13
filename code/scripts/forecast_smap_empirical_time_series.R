# For each point in time for a given empirical time series infer the
# Jacobian matrix with the S-map and perform abundance forecasts with
# an LSTM neural network

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/smap_jacobian.R")
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
source_python('code/functions/forecast_function.py')
source_python('code/functions/single_forecast_lstm.py')

# loading time series and defining settings ------------------------------ 
# proportion of data as training set
training_prop <- 0.7
# how many steps ahead to forecast
steps_ahead <- 3
# whether to normalize data for s-map
normalize <- FALSE
# whether to do cross-validation for LSTM hyperparameters
cv_LSTM <- FALSE
# empirical data (beninca_2009 or beninca_2015)
data <- "beninca_2015"
# kernel parameter for s-map
theta_seq <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 
               1, 2, 3, 4, 5, 6, 7, 8)
# load data
ts <- read.csv(paste("data/empirical_time_series/", data, ".csv", sep = ""), header = TRUE)
sp_names <- names(ts)[-1]
n_sp <- ncol(ts) - 1

# performing inference of Jacobian matrix and forecasts sequentially ------------------------------ 
# last training point to start with
start_last_train <- floor(nrow(ts) * training_prop)
# total number of forecasts
n_forecasts <- length((start_last_train + 1):nrow(ts)) - steps_ahead + 1
# objects to store results
smap_jacobians <- list()
smap_eigenvalue <- c()
smap_eigenvector <- list()
forecast_results <- list()
# for loop moving forward one point at a time
for (i in 1:n_forecasts) {
  print(i)
  # define current training and test sets
  training_set_curr <- i:(start_last_train - 1 + i)
  test_set_curr <- (start_last_train + i):(start_last_train + i + steps_ahead - 1)
  # prepare time series for s-map
  ts_smap <- ts[training_set_curr, ]
  names(ts_smap) <- c("time", paste("x", 1:n_sp, sep = ""))
  ts_smap_curr <- ts_smap
  # normalize training set for s-map
  if (normalize) {
    ts_smap_curr[ , -1] <- scale(ts_smap_curr[ , -1])
  }
  # select kernel parameter with cross-validation
  rmse_seq <- sapply(theta_seq, function(x) mean(smap_jacobian(x, ts = ts_smap_curr)[[2]]))
  theta <- theta_seq[which.min(rmse_seq)]
  # perform s-map with optimal theta
  smap_results <- smap_jacobian(ts_smap_curr, theta = theta)
  # s-map jacobians
  smap_J <- smap_results[[1]]
  # s-map intercepts
  smap_intercept <- smap_results[[3]]
  # jacobian at last time step
  smap_jacobians[[i]] <- smap_J[[length(smap_J)]]
  # eigendecomposition at last training point
  eigen_dec <- eigen(smap_J[[length(smap_J)]])
  # leading eigenvalue
  smap_eigenvalue[i] <- max(Re(eigen_dec$values))
  # leading eigenvector
  order_values <- order(Re(eigen_dec$values), decreasing = TRUE)
  smap_eigenvector[[i]] <- Re(eigen_dec$vectors)[ , order_values[1]]
  # perform forecasts
  forecast_results[[i]] <- single_forecast_lstm(training_data = as.matrix(ts_smap[ , -1]), 
                                                n_test = as.integer(steps_ahead), 
                                                cv = cv_LSTM)
}

# organize and save results ------------------------------ 
# merge s-map results into single data frame
time <- ts$time[start_last_train:(nrow(ts) - steps_ahead)]
smap_df <- as.data.frame(do.call("rbind", smap_eigenvector))
names(smap_df) <- paste("v1_", sp_names, sep = "")
lambda1 <- smap_eigenvalue
smap_df <- cbind(time, smap_df, lambda1)
# organize forecast data frame
last_train_time <- ts$time[rep(start_last_train:(nrow(ts) - steps_ahead), each = steps_ahead)]
time <- ts$time[unlist(sapply(seq(start_last_train + 1, nrow(ts) - steps_ahead + 1), 
                              seq, by = 1, length = steps_ahead, simplify = FALSE))]
forecast_df <- as.data.frame(do.call("rbind", forecast_results))
names(forecast_df) <- sp_names
forecast_df <- cbind(last_train_time, time, forecast_df)
# save r data
if (normalize) {
  save(smap_df, file = paste("results/forecast_analyses/sequential_smap_", data, "_",
                             "normalized_", training_prop, "_training_", steps_ahead, "_steps_ahead",
                             ".RData", sep = ""))
  save(smap_jacobians, file = paste("results/forecast_analyses/sequential_smap_jacobians_", data, "_",
                                    "normalized_", training_prop, "_training_", steps_ahead, "_steps_ahead",
                                    ".RData", sep = ""))
  save(forecast_df, file = paste("results/forecast_analyses/sequential_forecast_", data, "_",
                                 "normalized_", training_prop, "_training_", steps_ahead, "_steps_ahead",
                                 ".RData", sep = ""))
} else {
  save(smap_df, file = paste("results/forecast_analyses/sequential_smap_", data, "_",
                             training_prop, "_training_", steps_ahead, "_steps_ahead",
                             ".RData", sep = ""))
  save(smap_jacobians, file = paste("results/forecast_analyses/sequential_smap_jacobians_", data, "_",
                                    training_prop, "_training_", steps_ahead, "_steps_ahead",
                                    ".RData", sep = ""))
  save(forecast_df, file = paste("results/forecast_analyses/sequential_forecast_", data, "_",
                                 training_prop, "_training_", steps_ahead, "_steps_ahead",
                                 ".RData", sep = ""))
}
