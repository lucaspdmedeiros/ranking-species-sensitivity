# Code for Fig 2: application of the expected sensitivity ranking approach to 
# the 3-species food chain model 

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/food_chain.R")
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
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}

# plot settings ------------------------------
# whether to save plots
save_plots <- TRUE

# models and simulation settings ------------------------------
# number of species
n_sp <- 3
# model to use
func_name <- "food_chain"
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# amount of observational noise to add to time series
noise <- 0
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
# distribution of perturbed points 
pert_dist <- "gaussian"
# perturbation magnitude (percentage of average species standard deviation)
pert_magnitude <- 0.15
# time steps to evolve perturbed points
k <- 1
# whether k is fixed or variable
k_type <- "variable"

# prepare data to plot ------------------------------
# load results file
load(file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name, "_", n_sp, "_sp_",
                  pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                  k_type, "_time_step_k_", k, "_", jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
# extract expected and observed sensitivities
df_obs_exp_sens <- cbind(df[ , paste("s", 1:n_sp, sep = "")], 
                         df[ , paste("expected_sensitivity_", 1:n_sp, sep = "")])
# compute rank correlation between expected and observed sensitivities
cor_exp_sens <- apply(df_obs_exp_sens, 1, 
                      function(x) cor(as.numeric(x[1:n_sp]), 
                                      as.numeric(x[(n_sp+1):(2*n_sp)]), 
                                      method = "spearman"))
# percentage of points for each rank correlation value
table(cor_exp_sens) / length(cor_exp_sens)
# mean correlation value
mean(cor_exp_sens)
# add rank correlation to results data frame
df$spearman <- cor_exp_sens
# plot rank correlation over time
fig <- ggplot() +
  geom_line(data = df, aes(x = time, y = spearman), 
             size = 0.7, color = "gray80") +
  geom_point(data = df, aes(x = time, y = spearman), 
             size = 3.8, shape = 21, fill = "gray80") +
  geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed") +
  xlab("Time") +
  ylab(expression(atop("Correlation" ~ (rho) ~ "between expected",
                       "and observed sensitivities"))) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = c(0.11, 0.11))
# save plot
if (save_plots) {
  ggsave(paste("figs/fig2/fig2_", jacobian, "_J.pdf", sep = ""), 
         fig, width = 22, height = 10, units = "cm")
}
