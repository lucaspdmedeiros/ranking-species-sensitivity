# For each point along a time series apply random perturbations, evolve each perturbed point
# according to model equations and save results

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

# load synthetic time series ------------------------------
# number of species (2, 3, 4 or 5)
n_sp <- 2
# model to use (predator_prey, food_chain, lotka_volterra or consumer_resource)
func_name <- "predator_prey"
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# load result files
load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                  ts_length, "_points_", sampling_freq, "_sampling_freq_0_noise",
                  ".RData", sep = ""))
# load model settings
source("code/scripts/model_settings.R")
# if time series is stochastic, use deterministic model to evolve perturbed abundances
if (grepl("stochastic", func_name)) {
  stochastic_name <- func_name
  func_name <- substr(func_name, 1 , nchar(func_name)-11)
  source("code/scripts/model_settings.R")
  func_name <- stochastic_name
}

# perform perturbation simulations ------------------------------
# distribution of perturbed points (uniform, gaussian, or gaussian_proportional)
pert_dist <- "gaussian"
# perturbation magnitude (percentage of average species standard deviation)
pert_magnitude <- 0.15
scaled_pert_magnitude <- pert_magnitude * mean(apply(ts[ , -1], 2, sd))
# number of perturbations
n_pert <- 300
# time steps to evolve perturbed points (1 or 3)
k <- 1
# whether number of time points to evolve perturbations is fixed or variable
k_type <- "variable"
# time sequence to evolve perturbed points
times <- seq(0, k, by = time_step)
# data frames to store results
df_unpert <- data.frame()
df_full <- data.frame()
k_vec <- c()
# loop over all points in time series
for (j in 1:(nrow(ts)-1)) {
  print(j)
  state <- as.numeric(ts[j, -1])
  names(state) <- paste("x", 1:n_sp, sep = "")
  df_initial <- data.frame()
  df_final <- data.frame()
  # k inversely proportional to mean abundance rate of change
  if (k_type == "variable") {
    mean_rate_change <- mean(abs((as.numeric(ts[j+1, -1]) - as.numeric(ts[j, -1])) / as.numeric(ts[j, -1])))
    k_vec[j] <- round(k / sqrt(mean_rate_change), digits = 1)
    times_variable <- seq(0, k_vec[j], by = time_step)
    if (k_vec[j] == 0) {
      k_vec[j] <- 2 * time_step 
      times_variable <- seq(0, k_vec[j], by = time_step)
    }
    times_pert <- times_variable
  }
  # k remains the same for all points in time
  if (k_type == "fixed") {
    times_pert <- times
  }
  # apply n_pert perturbations (cloud of perturbations)
  for (i in 1:n_pert) {
    # sampling abundance perturbation
    if (pert_dist == "gaussian") {
      pert <- rnorm(n_sp, mean = 0, sd = scaled_pert_magnitude)
    }
    if (pert_dist == "gaussian_proportional") {
      pert <- rnorm(n_sp, mean = 0, sd = sqrt(state / sum(state)) * scaled_pert_magnitude)
    }
    if (pert_dist == "uniform") {
      x <- rnorm(n_sp, 0, 1)
      norm_x <- sqrt(sum(x^2))
      x_scaled <- (x / norm_x) * scaled_pert_magnitude
      u <- runif(1, 0, 1)^(1/n_sp)
      pert <- x_scaled * u
    }
    # applying perturbation
    state_initial <- state + pert
    # set perturbed abundance to zero, if negative
    state_initial[state_initial < 0] <- 0
    # add initial perturbed abundances to data frame
    df_initial <- rbind(df_initial, state_initial)
    # perform numerical integration
    ts_pert <- as.data.frame(ode(y = state_initial, times = times_pert, func = func, parms = parms, method = "ode45"))
    # add final perturbed abundances to data frame
    state_final <- as.numeric(ts_pert[nrow(ts_pert), -1])
    df_final <- rbind(df_final, state_final)
  }
  # obtain unperturbed final state
  ts_unpert <- as.data.frame(ode(y = state, times = times_pert, func = func, parms = parms, method = "ode45"))
  state_unpert <- as.numeric(ts_unpert[nrow(ts_unpert), -1])
  df_unpert <- rbind(df_unpert, state_unpert)
  # merge data frames
  names(df_initial) <- paste("x", 1:n_sp, sep = "")
  names(df_final) <- paste("x", 1:n_sp, sep = "")
  df_pert <- rbind(df_initial, df_final)
  df_pert$perturbation <- c(1:nrow(df_initial), 1:nrow(df_final))
  df_pert$type <- c(rep("initial", n_pert), rep("final", n_pert))
  df_pert$time <- j
  df_full <- rbind(df_full, df_pert)
}
# modifying data frames
ts <- ts[-nrow(ts), ]
df_unpert <- cbind(ts$time, df_unpert)
names(df_unpert) <- c("time", paste("x", 1:n_sp, sep = ""))

# saving objects as RData ------------------------------
if (k_type == "variable") {
  df_k <- data.frame(time = ts$time, k = k_vec)	
  save(df_k, file = paste("results/perturbation_analyses/time_step_k_", func_name, "_", n_sp, "_sp_",
                          k_type, "_time_step_k_", k, ".RData", sep = ""))
}
save(df_unpert, file = paste("results/perturbation_analyses/unperturbed_final_points_", func_name, "_", n_sp, "_sp_",
                             pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                             k_type, "_time_step_k_", k, ".RData", sep = ""))
save(df_full, file = paste("results/perturbation_analyses/perturbed_initial_final_points_", func_name, "_", n_sp, "_sp_", 
                           pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                           k_type, "_time_step_k_", k, ".RData", sep = ""))
