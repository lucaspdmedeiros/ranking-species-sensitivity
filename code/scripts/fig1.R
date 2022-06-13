# Code for Fig 1: illustration of time-varying species sensitivities to 
# perturbations for the 3-species food chain model

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
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}

# plot settings ------------------------------
# whether to save plots
save_plots <- TRUE
# set seed to replicate figure
set.seed(20)
# color palettes
palette1 <- c(pal_material("red")(9)[8], pal_material("amber")(9)[8])
palette2 <- brewer.pal(9, "Set1")[c(3, 2, 4)]

# simulate model to generate time series ------------------------------
# number of species
n_sp <- 3
# model to use
func_name <- "food_chain"
# final time series length
ts_length <- 3000
# amount of observational noise to add to time series
noise <- 0
# load model settings
source("code/scripts/model_settings.R")
# generate time series
ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
# remove transient
ts <- tail(ts, ts_length)
# simplifying times
ts$time <- seq(0, length = nrow(ts), by = time_step)

# perform perturbations ------------------------------
# perturbation magnitude (percentage of average species standard deviation)
pert_magnitude <- 0.15
scaled_pert_magnitude <- pert_magnitude * mean(apply(ts[ , -1], 2, sd))
# number of perturbations
n_pert <- 300
# time steps to evolve perturbed points
k <- 1
# time sequence to evolve perturbed points
times <- seq(0, k, by = time_step)
# define time points to apply perturbations
t_all <- c(145, 180)
# perform n_pert perturbations at the two different time points
df_full <- data.frame()
for (j in 1:length(t_all)) {
  # data frames to store results
  df_initial <- data.frame()
  df_final <- data.frame()
  # current abundances
  state <- as.numeric(ts[t_all[j], -1])
  names(state) <- paste("x", 1:n_sp, sep = "")
  # apply n_pert perturbations (cloud of perturbations)
  for (i in 1:n_pert) {
    # computing a gaussian perturbation around the current abundance
    pert <- rnorm(n_sp, mean = 0, sd = scaled_pert_magnitude)
    # perturbing abundances
    state_initial <- state + pert
    # set perturbed abundance to zero, if negative
    state_initial[state_initial < 0] <- 0
    # add initial perturbed abundances to data frame
    df_initial <- rbind(df_initial, state_initial)
    # evolve perturbed point over time
    ts_pert <- as.data.frame(ode(y = state_initial, times = times, func = func, parms = parms, method = "ode45"))
    # add final perturbed abundances to data frame
    state_final <- as.numeric(ts_pert[nrow(ts_pert), -1])
    df_final <- rbind(df_final, state_final)
  }
  # merge data frames
  names(df_initial) <- paste("x", 1:n_sp, sep = "")
  names(df_final) <- paste("x", 1:n_sp, sep = "")
  df_pert <- rbind(df_initial, df_final)
  df_pert$perturbation <- c(1:nrow(df_initial), 1:nrow(df_final))
  df_pert$type <- c(rep("initial", n_pert), rep("final", n_pert))
  df_pert$time <- t_all[j]
  df_full <- rbind(df_full, df_pert)
}

# plots for the two different time points ------------------------------
# define which species is most sensitive
sensitive <- c("sp2", "sp3")
for (i in 1:length(t_all)) {
  # defining time point
  t <- t_all[i]
  # plot 3D attractor with initial and final perturbation cloud
  fig <- plot_ly(x = ~x1, y = ~x2, z = ~x3, colors = palette1, showlegend = FALSE) %>% 
    add_trace(data = ts, type = 'scatter3d', mode = "lines",
              color = I('black'), line = list(width = 3)) %>% 
    add_markers(data = subset(df_full, time == t), color = ~type, 
                type = 'scatter3d', marker = list(size = 2.3)) %>% 
    layout(scene = list(xaxis = list(title = "Abundance species 1",
                                     titlefont = list(size = 22, 
                                                      family = "Arial, sans-serif",
                                                      color = palette2[1]), 
                                     tickfont = list(size = 14,
                                                     family = "Arial, sans-serif",
                                                     color = palette2[1]),
                                     range = c(0, 120),
                                     ticklen = 4.5,
                                     tickcolor = "#7F7F7F",
                                     gridcolor = "#7F7F7F",
                                     gridwidth = 0.6,
                                     zerolinewidth = 2.4),
                        yaxis = list(title = "Abundance species 2",
                                     titlefont = list(size = 22, 
                                                      family = "Arial, sans-serif",
                                                      color = palette2[2]),
                                     tickfont = list(size = 14,
                                                     family = "Arial, sans-serif",
                                                     color = palette2[2]),
                                     range = c(0, 120),
                                     ticklen = 4.5,
                                     tickcolor = "#7F7F7F",
                                     gridcolor = "#7F7F7F",
                                     gridwidth = 0.6,
                                     zerolinewidth = 2.4),
                        zaxis = list(title = "Abundance species 3",
                                     titlefont = list(size = 22, 
                                                      family = "Arial, sans-serif",
                                                      color = palette2[3]),
                                     tickfont = list(size = 14,
                                                     family = "Arial, sans-serif",
                                                     color = palette2[3]),
                                     range = c(0, 120),
                                     ticklen = 4.5,
                                     tickcolor = "#7F7F7F",
                                     gridcolor = "#7F7F7F",
                                     gridwidth = 0.6,
                                     zerolinewidth = 2.4),
                        camera = list(eye = list(x = 1.8, y = 0.8, z = 1)),
                        aspectratio = list(x = 0.8, y = 0.8, z = 0.8)),
           margin = list(l = 0,
                         r = 0,
                         b = 0,
                         t = 0))
  # save plot
  if (save_plots) {
    orca(fig, paste("figs/fig1/fig1_attractor_", sensitive[i], ".pdf", sep = ""), format = "pdf",
         width = 1000, height = 600)
  }
  # length of time series plot
  ts_plot_length <- 40
  # abundances at time t
  state <- as.numeric(ts[t, -1])
  names(state) <- paste("x", 1:n_sp, sep = "")
  # apply perturbation as an example
  state_initial <- state + c(7, 7, 7)
  # evolve perturbed point
  ts_pert <- as.data.frame(ode(y = state_initial, times = times, func = func, parms = parms, method = "ode45"))
  ts_pert$time <- seq(ts$time[t], length = nrow(ts_pert), by = time_step)
  # join both data frames with unperturbed and perturbed time series
  ts_all <- rbind(ts[(t-ts_plot_length):(t + nrow(ts_pert) - 1), ], ts_pert)
  ts_all$type <- "Unperturbed"
  ts_all$type[(nrow(ts_all) - nrow(ts_pert) + 1):nrow(ts_all)] <- "Perturbed"
  ts_all$type <- factor(ts_all$type, levels = c("Unperturbed", "Perturbed"))
  # building data frame for plot
  plot_df <- gather(ts_all, "species", "abundance", x1:x3)
  plot_df$species[plot_df$species == "x1"] <- "sp 1"
  plot_df$species[plot_df$species == "x2"] <- "sp 2"
  plot_df$species[plot_df$species == "x3"] <- "sp 3"
  # abundance time series plot with perturbed and unperturbed trajectories
  fig <- ggplot(data = plot_df, aes(x = time, y = abundance, color = species, linetype = type)) +
    geom_vline(xintercept = ts$time[t], size = 1.2, color = palette1[2]) +
    geom_vline(xintercept = ts$time[t+length(times) - 1], size = 1.2, color = palette1[1]) +
    geom_line(size = 2) +
    scale_color_manual(values = palette2) +
    xlab("Time") +
    ylab("Abundance") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size = 24),
          axis.text.x = element_text(size = 20),
          legend.position = "none")
  # save plot
  if (save_plots) {
    ggsave(paste("figs/fig1/fig1_time_series_", sensitive[i], ".pdf", sep = ""), 
           fig, width = 16, height = 8, units = "cm")
  }
}
