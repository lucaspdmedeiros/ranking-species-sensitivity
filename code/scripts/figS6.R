# Code for Fig S6: observed sensitivities, expected sensitivities, and
# eigenvector alignments for each species over time for each
# synthetic time series

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
source("code/functions/predator_prey.R")
source("code/functions/food_chain.R")
source("code/functions/consumer_resource.R")
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

# models and simulation settings ------------------------------
# number of species
n_sp <- c(2, 3, 3, 4, 5)
# model to use
func_name <- c("predator_prey", "food_chain", "lotka_volterra", 
               "lotka_volterra", "consumer_resource")
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# amount of observational noise to add to time series (0 or 0.1)
noise <- 0
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
# distribution of perturbed points (uniform, gaussian, or gaussian_proportional)
pert_dist <- "gaussian"
# perturbation magnitude (percentage of average species standard deviation)
pert_magnitude <- 0.15
# time steps to evolve perturbed points (1 or 3)
k <- 1
# whether k is fixed or variable
k_type <- "variable"
# whether k is misspecified or greatly_misspecified
k_error <- "none"
# whether to normalize data for s-map
normalize <- FALSE
# color palette
palette <- brewer.pal(9, "Set1")[c(3, 2, 4, 8, 9)]
# whether to save plots
save_plots <- TRUE

# plot observed sensitivity, expected sensitivity, and eigenvector alignment for each model ------------------------------
# loop over models
for (i in 1:length(func_name)) {
  # load results
  if (k_error == "none") {
    if (normalize) {
      load(file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name[i], "_", n_sp[i], "_sp_",
                        pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_normalized_",
                        jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
    } else {
      load(file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name[i], "_", n_sp[i], "_sp_",
                        pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_",
                        jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
    }
  }
  if (k_error == "misspecified") {
    load(file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name[i], "_", n_sp[i], "_sp_",
                      pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_misspecified_",
                      jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
  }
  if (k_error == "greatly_misspecified") {
    load(file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name[i], "_", n_sp[i], "_sp_",
                      pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_greatly_misspecified_",
                      jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
  }
  # data frames with (scaled) observed sensitivity
  df_obs_sens <- df[ , paste("s", 1:n_sp[i], sep = "")]
  df_obs_sens <- df_obs_sens / apply(df_obs_sens, 1, sum)
  df_obs_sens$time <- df$time
  names(df_obs_sens) <- c(paste("species", 1:n_sp[i], sep = " "), "time")
  df_obs_sens$type <- "Observed sensitivity"
  # data frames with (scaled) expected sensitivity
  df_expect_sens <- df[ , paste("expected_sensitivity_", 1:n_sp[i], sep = "")]
  df_expect_sens <- df_expect_sens / apply(df_expect_sens, 1, sum)
  df_expect_sens$time <- df$time
  names(df_expect_sens) <- c(paste("species", 1:n_sp[i], sep = " "), "time")
  df_expect_sens$type <- "Expected sensitivity"
  # data frames with (scaled) eigenvector alignment
  df_eigen_align <- df[ , paste("eigen_alignment_", 1:n_sp[i], sep = "")]
  df_eigen_align <- df_eigen_align / apply(df_eigen_align, 1, sum)
  df_eigen_align$time <- df$time
  names(df_eigen_align) <- c(paste("species", 1:n_sp[i], sep = " "), "time")
  df_eigen_align$type <- "Eigenvector alignment"
  # merge all data frames 
  full_df <- rbind(df_obs_sens, df_expect_sens, df_eigen_align)
  # creating data frame for plotting
  plot_df <- gather(data = full_df, key = "species", value = "variable", -c(time, type))
  plot_df$species <- as.factor(plot_df$species)
  plot_df$type <- factor(plot_df$type, levels = c("Observed sensitivity", 
                                                  "Expected sensitivity",
                                                  "Eigenvector alignment"))
  
  # plotting observed sensitivities and rankings over time ------------------------------ 
  # plot
  fig <- ggplot(data = plot_df, aes(x = time, y = variable, fill = species)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = palette) +
    facet_wrap(~type, ncol = 3) +
    xlab("Time") +
    ylab("Scaled variable") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          strip.text = element_text(size = 15),
          strip.background = element_rect(fill = "white", size = 1.5),
          axis.text.y = element_text(size = 13),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.key.size = unit(0.5, "cm"),
          legend.position = "top")
  # save plot
  if (save_plots) {
    ggsave(paste("figs/fig_sensitivity_rankings_over_time_", func_name[i], "_", n_sp[i], "_",
                 jacobian, "_", pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, 
                 "_time_step_k_", k, "_", noise, "_noise", ".pdf", sep = ""), 
           fig, width = 28, height = 8, units = "cm")
    
  }
}
