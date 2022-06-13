# Code for Fig S4: estimate leading Lyapunov vector for each point along 
# a time series and compute its alignment with the leading eigenvector

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
source("code/functions/predator_prey.R")
source("code/functions/food_chain.R")
source("code/functions/consumer_resource.R")
source("code/functions/original_and_tangent.R")
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(plotly)) {install.packages("plotly"); library(plotly)}
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}

# model and simulation settings ------------------------------
# to make results reproducible
set.seed(42)
# whether to save plots
save_plots <- TRUE
# number of species
n_sp_vec <- c(2, 3, 3, 4, 5)
# models to use
func_name_vec <- c("predator_prey", "food_chain", "lotka_volterra", 
                   "lotka_volterra", "consumer_resource")
models <- c("Predator-prey (2 sp)", "Food chain (3 sp)", "Food web (3 sp)",
            "Competitors (4 sp)", "Food web (5 sp)")
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# amount of observational noise added to time series
noise <- 0
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
# number of time steps to evolve perturbed points (1 or 3)
k <- 1
# whether number of time points to evolve perturbations is fixed or variable
k_type <- "variable"

# loop over models ------------------------------
plot_df <- data.frame()
for (j in 1:length(func_name_vec)) {
  print(models[j])
  # settings
  n_sp <- n_sp_vec[j]
  func_name <- func_name_vec[j]
  # load model settings
  source("code/scripts/model_settings.R")
  # load result files
  load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                    ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                    ".RData", sep = ""))
  if (k_type == "variable") {
    load(file = paste("results/perturbation_analyses/time_step_k_", func_name, "_", n_sp, "_sp_",
                      k_type, "_time_step_k_", k, ".RData", sep = ""))
  }
  ts <- ts[-nrow(ts), ]
  # analytical Jacobian matrix for each time step
  if (jacobian == "analytical") {
    J <- dlply(ts, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                     func = func,
                                                     parms = parms))
    # list with eigendecomposition of J matrices
    eigen_J_list <- lapply(J, eigen)
    # list to store leading eigenvector
    vectors <- list()
    for (i in 1:length(eigen_J_list)) {
      # decreasing order of eigenvalues
      values_order <- order(Re(eigen_J_list[[i]]$values), decreasing = TRUE)
      # extracting matrix with eigenvectors (same order as above)
      vectors[[i]] <- c(Re(eigen_J_list[[i]]$vectors)[ , values_order[1]])
    }
  }
  
  # loop over time series and compute leading Lyapunov vector ------------------------------
  # to store results
  cos_angle_eigenvector <- c()
  cos_angle_random <- c()
  lyap_exp_eigenvector <- c()
  lyap_exp_random <- c()
  # perturbation magnitude (percentage of average species standard deviation)
  scaled_pert_magnitude <- 0.05 * mean(apply(ts[ , -1], 2, sd))
  # loop over time series
  for (i in 1:nrow(ts)) {
    print(i)
    # number of integration time steps
    if (k_type == "variable") {
      k_curr <- df_k$k[i]
    } 
    if (k_type == "fixed") {
      k_curr <- 3
    }
    # time sequence
    times <- seq(0, k_curr, by = time_step)
    # current state
    state <- as.numeric(ts[i , -1])
    names(state) <- paste("x", 1:n_sp, sep = "")
    
    # tangent vector in the direction of leading eigenvector ------------------------------
    delta <- vectors[[i]]
    delta <- delta / sqrt(sum(delta^2))
    delta <- delta * scaled_pert_magnitude
    # fix potential negative perturbed abundances
    if (any(state + delta < 0)) {    
      delta <- -delta
      if (any(state + delta < 0)) {
        print(paste("t =", i))
        print("One or more perturbed abundances are negative")
      }
    }
    names(delta) <- paste("delta", 1:n_sp, sep = "")
    # merge the two vectors
    state_full <- c(state, delta)
    # integrate original and tangent dynamics simultaneously (normalizing tangent vector once)
    ts_original_tangent <- as.data.frame(ode(y = state_full, times = times, func = original_and_tangent, parms = parms, 
                                             original_func = func, n_sp = n_sp, method = "ode45"))
    # extract leading Lyapunov vector
    lyap_vector <- as.numeric(ts_original_tangent[nrow(ts_original_tangent), names(delta)])
    # extract leading Lyapunov exponent
    lyap_exp_eigenvector[i] <- (1/k_curr) * log(sqrt(sum(lyap_vector^2)) / scaled_pert_magnitude)
    # compute absolute value of cosine of angle between perturbation vector and and Lyapunov vector
    cos_angle_eigenvector[i] <- abs(sum(delta * lyap_vector) / (sqrt(sum(delta^2)) * sqrt(sum(lyap_vector^2))))
    
    # tangent vector in random direction ------------------------------
    delta <- rnorm(n_sp, 0, 1)
    delta <- delta / sqrt(sum(delta^2))
    delta <- delta * scaled_pert_magnitude
    # fix potential negative perturbed abundances
    if (any(state + delta < 0)) {
      while(any(state + delta < 0)) {
        delta <- rnorm(n_sp, 0, 1)
        delta <- delta / sqrt(sum(delta^2))
        delta <- delta * scaled_pert_magnitude
      }
    }
    names(delta) <- paste("delta", 1:n_sp, sep = "")
    # merge the two vectors
    state_full <- c(state, delta)
    # integrate original and tangent dynamics simultaneously (normalizing tangent vector once)
    ts_original_tangent <- as.data.frame(ode(y = state_full, times = times, func = original_and_tangent, parms = parms, 
                                             original_func = func, n_sp = n_sp, method = "ode45"))
    # extract leading Lyapunov vector
    lyap_vector <- as.numeric(ts_original_tangent[nrow(ts_original_tangent), names(delta)])
    # extract leading Lyapunov exponent
    lyap_exp_random[i] <- (1/k_curr) * log(sqrt(sum(lyap_vector^2)) / scaled_pert_magnitude)
    # compute absolute value of cosine of angle between perturbation vector and and Lyapunov vector
    cos_angle_random[i] <- abs(sum(delta * lyap_vector) / (sqrt(sum(delta^2)) * sqrt(sum(lyap_vector^2))))
  }
  
  # build data frame with results ------------------------------
  results_df <- data.frame(time = rep(ts$time, 2),
                           cos_angle = c(cos_angle_eigenvector, cos_angle_random),
                           lyap_exp = c(lyap_exp_eigenvector, lyap_exp_random),
                           type = c(rep("Eigenvector", length(cos_angle_eigenvector)),
                                    rep("Random", length(cos_angle_random))),
                           model = models[j])
  # merge data frames
  plot_df <- rbind(plot_df, results_df)
}

# plotting distribution of angle cosine for eigenvector and random perturbations ------------------------------
# changing variables to factor
plot_df$model <- factor(plot_df$model, levels = c("Predator-prey (2 sp)", "Food chain (3 sp)", "Food web (3 sp)",
                                                  "Competitors (4 sp)", "Food web (5 sp)"))
plot_df$type <- factor(plot_df$type, levels = c("Eigenvector", "Random"))
# mean values of angle cosines
ddply(plot_df, c("model", "type"), summarise,
      mean_cos_angle = mean(cos_angle))
# plot
fig <- ggplot(data = plot_df, aes(x = type, y = cos_angle, group = type)) +
  geom_boxplot(size = 1.2, outlier.size = 1.5) +
  facet_wrap(~model, nrow = 3) +
  xlab(TeX("Direction of initial perturbation ($\\mathbf{p}(t)$)")) +
  ylab(TeX("Cosine of the angle between $\\mathbf{p}(t)$ and $\\mathbf{p}(t+k)$")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 21),
        axis.text.x = element_text(size = 16),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/fig_cosine_v1_lyapunov_vector_", k_type, "_time_step_k_", 
               k, ".pdf", sep = ""), 
         fig, width = 20, height = 30, units = "cm")
}
