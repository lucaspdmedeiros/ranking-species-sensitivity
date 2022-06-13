# For each point in time for a given synthetic time series compute:
# (1) Jacobian matrix, (2) eigenvector alignments, 
# (3) expected sensitivities, and (4) observed sensitivities

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
source("code/functions/predator_prey.R")
source("code/functions/food_chain.R")
source("code/functions/consumer_resource.R")
source("code/functions/lotka_volterra_stochastic.R")
source("code/functions/predator_prey_stochastic.R")
source("code/functions/food_chain_stochastic.R")
source("code/functions/consumer_resource_stochastic.R")
source("code/functions/smap_jacobian.R")
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
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
if(!require(expm)) {install.packages("expm"); library(expm)}
if(!require(glmnet)) {install.packages("glmnet"); library(glmnet)}

# model and simulation settings ------------------------------
# number of species (2, 3, 4 or 5)
n_sp <- 2
# model to use (predator_prey, food_chain, lotka_volterra or consumer_resource)
func_name <- "predator_prey"
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# amount of observational noise added to time series (0 or 0.1)
noise <- 0
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
# whether to normalize data for s-map
normalize <- FALSE
# fraction of time series to train s-map
frac_train <- 0.5
# kernel parameter for s-map
theta_seq <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 
               1, 2, 3, 4, 5, 6, 7, 8)
# load model settings
source("code/scripts/model_settings.R")
# distribution of perturbed points (uniform, gaussian, or gaussian_proportional)
pert_dist <- "gaussian"
# perturbation magnitude (percentage of average species standard deviation)
pert_magnitude <- 0.15
# number of time steps to evolve perturbed points (1 or 3)
k <- 1
# whether number of time points to evolve perturbations is fixed or variable
k_type <- "variable"
# whether k is misspecified or greatly_misspecified (default is none)
k_error <- "none"
# load result files
load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                  ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                  ".RData", sep = ""))
load(file = paste("results/perturbation_analyses/unperturbed_final_points_", func_name, "_", n_sp, "_sp_",
                  pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                  k_type, "_time_step_k_", k, ".RData", sep = ""))
load(file = paste("results/perturbation_analyses/perturbed_initial_final_points_", func_name, "_", n_sp, "_sp_", 
                  pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                  k_type, "_time_step_k_", k, ".RData", sep = ""))
if (k_type == "variable") {
  load(file = paste("results/perturbation_analyses/time_step_k_", func_name, "_", n_sp, "_sp_",
                    k_type, "_time_step_k_", k, ".RData", sep = ""))
}
# removing last time series point that is not present in results data frames
ts <- ts[-nrow(ts), ]

# compute Jacobian matrix ------------------------------ 
# list to store results
J <- list()
# time points to train s-map
train_points <- 1:ceiling(nrow(ts) * frac_train)
# time points to infer Jacobians
test_points <- ceiling(nrow(ts) * frac_train):nrow(ts)
# time series to train s-map
ts_train <- ts[train_points, ]
# time series to infer Jacobians
ts <- ts[test_points, ]
# subset other data frames
df_full <- df_full[!is.na(match(df_full$time, test_points)), ]
df_unpert <- df_unpert[!is.na(match(df_unpert$time, test_points)), ]
# analytical Jacobian matrix for each time
if (jacobian == "analytical") {
  J <- dlply(ts, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                   func = func,
                                                   parms = parms))
}
# infer Jacobian matrix sequentially with the s-map using only past data
if (jacobian == "smap") {
  # infer Jacobian sequentially
  for (i in 1:nrow(ts)) {
    print(i)
    ts_train_curr <- ts_train
    # normalize time series
    if (normalize) {
      ts_train_curr[ , -1] <- scale(ts_train_curr[ , -1])
    }
    # select kernel parameter using cross validation
    rmse_seq <- sapply(theta_seq, function(x) mean(smap_jacobian(x, ts = ts_train_curr)[[2]]))
    theta <- theta_seq[which.min(rmse_seq)]
    # fit smap to past points
    J_train <- smap_jacobian(ts_train_curr, theta = theta)[[1]]
    # save last jacobian
    J[[i]] <- J_train[[length(J_train)]]
    # add a point to training data
    ts_train <- rbind(ts_train, ts[i + 1, ])
    # remove first point of training data
    ts_train <- ts_train[-1, ]
  }
}

# compute eigenvector alignments and expected sensitivities ------------------------------ 
# list with eigendecomposition of J matrices
eigen_J_list <- lapply(J, eigen)
# objects to save results
values <- list()
vectors <- list()
eigen_alignments <- list()
expected_sensitivities <- list()
for (i in 1:length(J)) {
  print(i)
  # current time point
  t <- ts$time[i]
  # decreasing order of eigenvalues
  values_order <- order(Re(eigen_J_list[[i]]$values), decreasing = TRUE)
  # extracting ordered eigenvalues
  values[[i]] <- Re(eigen_J_list[[i]]$values)[values_order]
  # extracting leading eigenvector 
  vectors[[i]] <- Re(eigen_J_list[[i]]$vectors)[ , values_order[1]]
  # eigenvector alignments
  eigen_alignments[[i]] <- abs(vectors[[i]] / sqrt(sum(vectors[[i]]^2)))
  # standard deviation of perturbations (used for expected sensitivities)
  pert_sd <- rep(1, n_sp)
  # value of k (time step for which perturbations were evolved)
  if (k_type == "fixed") {
    k_curr <- k
  }
  if (k_type == "variable") {
    k_curr <- df_k$k[t]
    if (k_error == "misspecified") {
      k_curr <- k
    }
  }
  if (k_error == "greatly_misspecified") {
    k_curr <- k_curr + rnorm(1, 0, k_curr)
    if (k_curr < 0) {
      k_curr <- 0.1
    }
    pert_sd <- rep(1, n_sp) + rnorm(n_sp, 0, 1)
    pert_sd[pert_sd < 0] <- 0.001
  }
  # initial covariance matrix of perturbations
  Sigma_initial <- diag(pert_sd^2)
  # exponential of k*J
  M <- expm(k_curr * J[[i]])
  # final covariance matrix of perturbations
  Sigma_final <- M %*% Sigma_initial %*% t(M) 
  # expected sensitivities
  expected_sensitivities[[i]] <- diag(Sigma_final)
}

# organize results into single data frame ------------------------------ 
# data frame with eigenvalues
df_values <- data.frame(matrix(unlist(values), nrow = length(values), byrow = TRUE))
names(df_values) <- paste("lambda", 1:n_sp, sep = "_")
# data frame with eigenvector alignments
df_eigen_alignments <- data.frame(matrix(unlist(eigen_alignments), nrow = length(eigen_alignments), byrow = TRUE))
names(df_eigen_alignments) <- paste("eigen_alignment", 1:n_sp, sep = "_")
# data frame with expected sensitivities 
df_expected_sensitivities <- abs(data.frame(matrix(unlist(expected_sensitivities), 
                                                   nrow = length(expected_sensitivities), byrow = TRUE)))
names(df_expected_sensitivities) <- paste("expected_sensitivity", 1:n_sp, sep = "_")
# data frame with time
df_time <- data.frame(time = ts$time)
# data frame with time series
df_abund <- ts[ , -1]
# data frame with rate of change
percent_change <- abs(rbind((ts[2:nrow(ts), 2:(n_sp+1)] - ts[1:(nrow(ts)-1), 2:(n_sp+1)]) / 
                              ts[1:(nrow(ts)-1), 2:(n_sp+1)]))
df_percent_change <- rbind(rep(NA, n_sp), percent_change)
names(df_percent_change) <- paste("delta", paste("x", 1:n_sp, sep = ""), sep = "_")
# merging all results into a single data frame
df <- cbind(df_time, df_abund, df_values, df_eigen_alignments, df_expected_sensitivities, df_percent_change)

# compute average species sensitivities to perturbations ------------------------------ 
s1 <- c()
s2 <- c()
s3 <- c()
s4 <- c()
s5 <- c()
df_initial <- df_full[df_full$type == "initial", ]
df_final <- df_full[df_full$type == "final", ]
# if noise is not zero, load time series without noise
if (noise != 0) {
  load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                    ts_length, "_points_", sampling_freq, "_sampling_freq_0_noise",
                    ".RData", sep = ""))
  ts <- ts[test_points, ]
}
# compute sensitivities
for (i in 1:nrow(df)) {
  t <- ts$time[i]
  s1[i] <- mean((df_final[df_final$time == t, "x1"] - df_unpert[i, "x1"])^2) / 
    mean((df_initial[df_initial$time == t, "x1"] - ts[i, "x1"])^2)
  s2[i] <- mean((df_final[df_final$time == t, "x2"] - df_unpert[i, "x2"])^2) / 
    mean((df_initial[df_initial$time == t, "x2"] - ts[i, "x2"])^2)
  if (n_sp == 3) {
    s3[i] <- mean((df_final[df_final$time == t, "x3"] - df_unpert[i, "x3"])^2) / 
      mean((df_initial[df_initial$time == t, "x3"] - ts[i, "x3"])^2)
  }
  if (n_sp == 4) {
    s3[i] <- mean((df_final[df_final$time == t, "x3"] - df_unpert[i, "x3"])^2) / 
      mean((df_initial[df_initial$time == t, "x3"] - ts[i, "x3"])^2)
    s4[i] <- mean((df_final[df_final$time == t, "x4"] - df_unpert[i, "x4"])^2) / 
      mean((df_initial[df_initial$time == t, "x4"] - ts[i, "x4"])^2)
  }
  if (n_sp == 5) {
    s3[i] <- mean((df_final[df_final$time == t, "x3"] - df_unpert[i, "x3"])^2) / 
      mean((df_initial[df_initial$time == t, "x3"] - ts[i, "x3"])^2)
    s4[i] <- mean((df_final[df_final$time == t, "x4"] - df_unpert[i, "x4"])^2) / 
      mean((df_initial[df_initial$time == t, "x4"] - ts[i, "x4"])^2)
    s5[i] <- mean((df_final[df_final$time == t, "x5"] - df_unpert[i, "x5"])^2) / 
      mean((df_initial[df_initial$time == t, "x5"] - ts[i, "x5"])^2)
  }
}
df$s1 <- s1
df$s2 <- s2
df$s3 <- s3
df$s4 <- s4
df$s5 <- s5

# saving data frame ------------------------------ 
# save results
if (k_error == "none") {
  if (normalize) {
    save(df, file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name, "_", n_sp, "_sp_",
                          pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_normalized_",
                          jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
  } else {
    save(df, file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name, "_", n_sp, "_sp_",
                          pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_",
                          jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
  }
}
if (k_error == "misspecified") {
  save(df, file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name, "_", n_sp, "_sp_",
                        pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_misspecified_",
                        jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
}
if (k_error == "greatly_misspecified") {
  save(df, file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name, "_", n_sp, "_sp_",
                        pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_greatly_misspecified_",
                        jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
}
