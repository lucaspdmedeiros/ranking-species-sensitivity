# Code for Fig S1-S3: illustration of the expected sensitivity and eigenvector
# ranking approaches under Lotka-Volterra dynamics at equilibrium

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
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
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}
if(!require(expm)) {install.packages("expm"); library(expm)}

# models and simulation settings ------------------------------
# to make results reproducible
set.seed(4)
# whether to save plots
save_plots <- TRUE
# number of species
n_sp <- 3
# model to use
func <- lotka_volterra
# which case to use
case <- "case_1"
# LV parameters for lambda_1 > 0, lambda_2 < 0, and lambda_3 < 0
if (case == "case_1") {
  A <- matrix(c(1, -2, 0, 
                0, -1, 0,
                0, 2, -3), nrow = 3, ncol = 3, byrow = TRUE)
  r <- c(-A %*% c(1, 1, 1))
}
# LV parameters for lambda_1 > 0, lambda_2 > 0, and lambda_3 < 0
if (case == "case_2") {
  A <- matrix(c(4, 1/2, 0, 
                1/2, -10, -8,
                0, -8, 1), nrow = 3, ncol = 3, byrow = TRUE)
  r <- c(-A %*% c(1, 1, 1))
}
# LV parameters for lambda_1 > 0 (complex), lambda_2 > 0 (complex), and lambda_3 < 0
if (case == "case_3") {
  A <- matrix(c(-4, -3, 2,
                -2, 1, 2,
                5, 2, 0), nrow = 3, ncol = 3, byrow = TRUE)
  r <- c(-A %*% c(1, 1, 1))
}
# merge parameters
parms <- c(r, A)
names(parms) <- c(paste("r", 1:n_sp, sep = ""),
                  paste(paste("a", 1:n_sp, sep = ""), 
                        rep(1:n_sp, each = n_sp), sep = ""))

# eigenvalues and eigenvectors ------------------------------
# equilibrium abundances
N <- solve(A, -r)
# Jacobian at equilibrium
J <- diag(N) %*% A
# eigendecomposition
eigendec <- eigen(J)
# ordered eigenvalues
order_values <- order(Re(eigendec$values), decreasing = TRUE)
eigenvalues <- eigendec$values[order_values]
print(eigenvalues)
# ordered eigenvectors
eigenvectors <- eigendec$vectors[ , order_values]
print(eigenvectors[ , 1])

# expected sensitivities ------------------------------
# time steps to use
k <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
# standard deviation of Gaussian perturbation distribution
pert_sd <- rep(0.05 * mean(N), n_sp)
# initial covariance matrix of perturbations
Sigma_initial <- diag(pert_sd^2)
# compute expected sensitivities
M <- list()
expected_sensitivities <- list()
for (i in 1:(length(k)-1)) {
  # exponential of k*J
  M[[i]] <- expm(k[i+1] * J)
  # final covariance matrix of perturbations
  Sigma_final <- M[[i]] %*% Sigma_initial %*% t(M[[i]]) 
  # expected sensitivities
  expected_sensitivities[[i]] <- diag(Sigma_final)
  print(order(expected_sensitivities[[i]], decreasing = TRUE))
}

# perturbation analyses ------------------------------
# perform gaussian perturbations
n_points <- 2000
perturbed_abund <- list()
for (i in 1:n_points) {
  pert <- rnorm(n = n_sp, mean = 0, sd = 0.05 * mean(N))
  perturbed_abund[[i]] <- N + pert
}
# to integrate perturbed abundances
time_step <- 0.01
time_final <- 0.5
times <- seq(0, time_final, by = time_step)
# integrate perturbed abundances
results_df <- data.frame()
for (i in 1:length(perturbed_abund)) {
  ts <- as.data.frame(ode(y = perturbed_abund[[i]], times = times, func = func, parms = parms, method = "ode45"))
  final_abund <- ts[!is.na(match(ts$time, k)), ]
  names(final_abund) <- c("time", paste("x", 1:n_sp, sep = ""))
  final_abund$perturbation <- i
  results_df <- rbind(results_df, final_abund)
}

# plotting cloud of perturbations over time ------------------------------
# create plot data frame
plot_df <- subset(results_df, time == 0.4)
plot_df$k <- as.factor(plot_df$time)
min_scale <- min(plot_df[ , c("x1", "x2", "x3")])
max_scale <- max(plot_df[ , c("x1", "x2", "x3")])
# plot for sp 1 and sp 2
fig <- ggplot(data = plot_df, aes(x = x1, y = x2)) +
  geom_point(size = 2, alpha = 0.5, color = "gray60") +
  geom_point(data = data.frame(x1 = N[1], x2 = N[2], x3 = N[3]), 
             aes(x = x1, y = x2), size = 4, color = "black") +
  xlab(TeX("Abundance sp 1 ($N_1$)")) +
  ylab(TeX("Abundance sp 2 ($N_2$)")) +
  scale_x_continuous(limits = c(min_scale, max_scale)) +
  scale_y_continuous(limits = c(min_scale, max_scale)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 15),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/fig_perturbation_evolution_equilibrium_LV_sp1_sp2_",
               case, ".pdf", sep = ""), 
         fig, width = 10, height = 10, units = "cm")
}
# plot for sp 1 and sp 3
fig <- ggplot(data = plot_df, aes(x = x1, y = x3)) +
  geom_point(size = 2, alpha = 0.5, color = "gray60") +
  geom_point(data = data.frame(x1 = N[1], x2 = N[2], x3 = N[3]), 
             aes(x = x1, y = x3), size = 4, color = "black") +
  xlab(TeX("Abundance sp 1 ($N_1$)")) +
  ylab(TeX("Abundance sp 3 ($N_3$)")) +
  scale_x_continuous(limits = c(min_scale, max_scale)) +
  scale_y_continuous(limits = c(min_scale, max_scale)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 15),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/fig_perturbation_evolution_equilibrium_LV_sp1_sp3_",
               case, ".pdf", sep = ""), 
         fig, width = 10, height = 10, units = "cm")
}
# plot for sp 2 and sp 3
fig <- ggplot(data = plot_df, aes(x = x2, y = x3)) +
  geom_point(size = 2, alpha = 0.5, color = "gray60") +
  geom_point(data = data.frame(x1 = N[1], x2 = N[2], x3 = N[3]), 
             aes(x = x2, y = x3), size = 4, color = "black") +
  xlab(TeX("Abundance sp 2 ($N_2$)")) +
  ylab(TeX("Abundance sp 3 ($N_3$)")) +
  scale_x_continuous(limits = c(min_scale, max_scale)) +
  scale_y_continuous(limits = c(min_scale, max_scale)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 15),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/fig_perturbation_evolution_equilibrium_LV_sp2_sp3_",
               case, ".pdf", sep = ""), 
         fig, width = 10, height = 10, units = "cm")
}

# plotting species sensitivities for different time steps k ------------------------------
# function to compute species sensitivity (mean squared deviation from equilibrium point)
msd <- function(x, avg) {
  mean((x - avg)^2, na.rm = TRUE)
}
# computing species sensitivity for each value of k
summ_df <- ddply(results_df, c("time"), summarise,
                 x1 = msd(x1, avg = N[1]),
                 x2 = msd(x2, avg = N[2]),
                 x3 = msd(x3, avg = N[3]))
summ_df <- gather(summ_df, "species", "sensitivity", -time)
summ_df$species[summ_df$species == "x1"] <- "Species 1"
summ_df$species[summ_df$species == "x2"] <- "Species 2"
summ_df$species[summ_df$species == "x3"] <- "Species 3"
# color palette
palette <- brewer.pal(9, "Set1")[c(3, 2, 4)]
# plot of sensitivity of each species as a function of k
fig <- ggplot(data = summ_df, aes(x = time, y = sensitivity, color = species)) +
  geom_point(size = 4) +
  scale_color_manual(values = palette) +
  xlab("Time") +
  ylab("Species sensitivity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 17),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.8, "cm"),
        legend.position = "bottom")
if (save_plots) {
  ggsave(paste("figs/fig_sensitivity_vs_k_equilibrium_LV_",
               case, ".pdf", sep = ""), 
         fig, width = 14, height = 14, units = "cm")
}
