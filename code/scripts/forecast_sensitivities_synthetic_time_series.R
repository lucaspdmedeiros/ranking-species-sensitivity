# For each point in time for a given synthetic time series compute:
# (1) Jacobian matrix, (2) eigenvector alignments, 
# (3) expected sensitivities, and (4) average forecast errors

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
noise <- 0.1
# whether to use analytical or smap Jacobian
jacobian <- "smap"
# whether to normalize data for s-map
normalize <- FALSE
# fraction of time series to train s-map
frac_train <- 0.5
# how many steps ahead to forecast
steps_ahead <- 3
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
k <- steps_ahead
# whether number of time points to evolve perturbations is fixed or variable
k_type <- "fixed"
# load time series and perturbations analysis results
load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                  ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                  ".RData", sep = ""))
load(file = paste("results/perturbation_analyses/unperturbed_final_points_", func_name, "_", n_sp, "_sp_",
                  pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                  k_type, "_time_step_k_", k, ".RData", sep = ""))
load(file = paste("results/perturbation_analyses/perturbed_initial_final_points_", func_name, "_", n_sp, "_sp_", 
                  pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                  k_type, "_time_step_k_", k, ".RData", sep = ""))
# load s-map and forecast results
if (normalize) {
  load(paste("results/forecast_analyses/sequential_smap_jacobians_", func_name, "_", n_sp,
             "_normalized_", frac_train, "_training_", steps_ahead, "_steps_ahead",
             ".RData", sep = ""))
  load(paste("results/forecast_analyses/sequential_forecast_", func_name, "_", n_sp, 
             "_normalized_", frac_train, "_training_", steps_ahead, "_steps_ahead",
             ".RData", sep = ""))
} else {
  load(paste("results/forecast_analyses/sequential_smap_jacobians_", func_name, "_", n_sp, "_",
             frac_train, "_training_", steps_ahead, "_steps_ahead",
             ".RData", sep = ""))
  load(paste("results/forecast_analyses/sequential_forecast_", func_name, "_", n_sp, "_",
             frac_train, "_training_", steps_ahead, "_steps_ahead",
             ".RData", sep = ""))
}

# compute eigenvector alignments and expected sensitivities ------------------------------ 
# extract test set time points used for forecasting
forecast_steps_ahead <- forecast_df[seq(steps_ahead, nrow(forecast_df), by = steps_ahead), ]
test_points <- forecast_steps_ahead$last_train_time
# subset time series
ts <- ts[test_points, ]
# subset other data frames
df_full <- df_full[!is.na(match(df_full$time, test_points)), ]
df_unpert <- df_unpert[!is.na(match(df_unpert$time, test_points)), ]
# list with eigendecomposition of J matrices
J <- smap_jacobians
eigen_J_list <- lapply(J, eigen)
# objects to save results
values <- list()
vectors <- list()
eigen_alignments <- list()
expected_sensitivities <- list()
for (i in 1:length(J)) {
  print(i)
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
  k_curr <- steps_ahead
  # initial covariance matrix of perturbations
  Sigma_initial <- diag(pert_sd)
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

# compute average forecast errors ------------------------------ 
epsilon_1 <- c()
epsilon_2 <- c()
epsilon_3 <- c()
epsilon_4 <- c()
epsilon_5 <- c()
df_initial <- df_full[df_full$type == "initial", ]
df_final <- df_full[df_full$type == "final", ]
# compute average forecast errors
for (i in 1:nrow(forecast_steps_ahead)) {
  t <- forecast_steps_ahead$last_train_time[i]
  epsilon_1[i] <- mean(sqrt((df_final[df_final$time == t, "x1"] - forecast_steps_ahead[i, "x1"])^2))
  epsilon_2[i] <- mean(sqrt((df_final[df_final$time == t, "x2"] - forecast_steps_ahead[i, "x2"])^2))
  if (n_sp == 3) {
    epsilon_3[i] <- mean(sqrt((df_final[df_final$time == t, "x3"] - forecast_steps_ahead[i, "x3"])^2))
  }
  if (n_sp == 4) {
    epsilon_3[i] <- mean(sqrt((df_final[df_final$time == t, "x3"] - forecast_steps_ahead[i, "x3"])^2))
    epsilon_4[i] <- mean(sqrt((df_final[df_final$time == t, "x4"] - forecast_steps_ahead[i, "x4"])^2))
  }
  if (n_sp == 5) {
    epsilon_3[i] <- mean(sqrt((df_final[df_final$time == t, "x3"] - forecast_steps_ahead[i, "x3"])^2))
    epsilon_4[i] <- mean(sqrt((df_final[df_final$time == t, "x4"] - forecast_steps_ahead[i, "x4"])^2))
    epsilon_5[i] <- mean(sqrt((df_final[df_final$time == t, "x5"] - forecast_steps_ahead[i, "x5"])^2))
  }
}
df$epsilon_1 <- epsilon_1
df$epsilon_2 <- epsilon_2
df$epsilon_3 <- epsilon_3
df$epsilon_4 <- epsilon_4
df$epsilon_5 <- epsilon_5

# saving data frame ------------------------------ 
# save results
if (normalize) {
  save(df, file = paste("results/forecast_analyses/forecasts_sensitivities_", func_name, "_", n_sp, "_sp_",
                        pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_normalized_",
                        jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
} else {
  save(df, file = paste("results/forecast_analyses/forecasts_sensitivities_", func_name, "_", n_sp, "_sp_",
                        pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_",
                        jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
}
