# Code for figure 4A (also Figs S21-S25): 
# species expected sensitivities over time for each empirical time series

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

# plot settings ------------------------------
# whether to save plots
save_plots <- TRUE
# color palettes
palette <- brewer.pal(9, "Set1")[c(3, 2, 4, 8)]
# proportion of data as training set
training_prop <- 0.7
# how many steps ahead forecasted
steps_ahead <- 3
# empirical data sets to use
data <- c("beninca_2009", "beninca_2015")
# time steps of perturbed points
k <- steps_ahead
# whether to normalize data for s-map
normalize <- FALSE

# plotting time series colored by species sensitivity over time ------------------------------ 
# data frame to store results
full_results_df <- data.frame()
for (i in 1:length(data)) {
  # load empirical time series
  ts <- read.csv(paste("data/empirical_time_series/", data[i], ".csv", sep = ""), header = TRUE)
  sp_names <- names(ts)[-1]
  n_sp <- ncol(ts) - 1
  # load s-map results
  if (normalize) {
    load(paste("results/forecast_analyses/sequential_smap_", data[i], "_normalized_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
    load(paste("results/forecast_analyses/sequential_smap_jacobians_", data[i], "_normalized_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
  } else {
    load(paste("results/forecast_analyses/sequential_smap_", data[i], "_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
    load(paste("results/forecast_analyses/sequential_smap_jacobians_", data[i], "_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
  }
  
  # computing expected species sensitivities ------------------------------ 
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
  
  # plots ------------------------------ 
  # last training point to start with
  start_last_train <- which(ts$time == smap_df$time[1])
  # first test set point
  start_first_test <- which(ts$time == smap_df$time[2])
  # building data frame for plotting
  exp_sens_plot_df <- rbind(matrix(NA, nrow = start_last_train - 1, ncol = n_sp),
                            matrix(unlist(expected_sensitivities), nrow = length(expected_sensitivities), 
                                   byrow = TRUE),
                            matrix(NA, nrow = steps_ahead, ncol = n_sp))
  colnames(exp_sens_plot_df) <- paste("exp_sens", names(ts)[-1], sep = "_")
  # maximum and minimum expected sensitivity values
  max_exp_sen <- max(exp_sens_plot_df, na.rm = TRUE)
  min_exp_sen <- min(exp_sens_plot_df, na.rm = TRUE)
  # data frame for plotting
  ts_plot <- cbind(ts, exp_sens_plot_df)
  # removing steps_ahead last points
  ts_plot <- ts_plot[1:(nrow(ts) - steps_ahead), ]
  # plots for both time series
  if (data[i] == "beninca_2009") {
    # adding days
    ts_plot$day <- ts_plot$time
    # time series plot for calanoids
    fig_calanoids <- ggplot(ts_plot, aes(x = day, y = calanoids, fill = exp_sens_calanoids)) +
      geom_ribbon(aes(xmin = day[1], xmax = day[start_last_train]),
                  color = "gray90", fill = "gray90") +
      geom_line(size = 0.5) + 
      geom_point(size = 3, shape = 21) + 
      scale_fill_gradient(low = "#FFFFFF", high = "#4DAF4A", 
                          limits = c(min_exp_sen, max_exp_sen)) +
      scale_y_continuous(limits = c(0, 2.9)) +
      xlab("") +
      ylab("") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.45, "cm"))
    if (save_plots) {
      ggsave(paste("figs/fig4/fig4_ts_beninca_2009_",
                   training_prop, "_training_", steps_ahead, "_steps_ahead",
                   "_calanoids.pdf", sep = ""), 
             fig_calanoids, width = 34, height = 5, units = "cm")
    }
    # time series plot for nanophytoplankton
    fig_nanophytoplankton <- ggplot(ts_plot, aes(x = day, y = nanophytoplankton, fill = exp_sens_nanophytoplankton)) +
      geom_ribbon(aes(xmin = day[1], xmax = day[start_last_train]),
                  color = "gray90", fill = "gray90") +
      geom_line(size = 0.5) + 
      geom_point(size = 3, shape = 21) + 
      scale_fill_gradient(low = "#FFFFFF", high = "#377EB8", 
                          limits = c(min_exp_sen, max_exp_sen)) +
      scale_y_continuous(limits = c(0, 2.9)) +
      xlab("") +
      ylab("") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.45, "cm"))
    if (save_plots) {
      ggsave(paste("figs/fig4/fig4_ts_beninca_2009_",
                   training_prop, "_training_", steps_ahead, "_steps_ahead",
                   "_nanophytoplankton.pdf", sep = ""), 
             fig_nanophytoplankton, width = 34, height = 5, units = "cm")
    }
    # time series plot for picophytoplankton
    fig_picophytoplankton <- ggplot(ts_plot, aes(x = day, y = picophytoplankton, fill = exp_sens_picophytoplankton)) +
      geom_ribbon(aes(xmin = day[1], xmax = day[start_last_train]),
                  color = "gray90", fill = "gray90") +
      geom_line(size = 0.5) + 
      geom_point(size = 3, shape = 21) + 
      scale_fill_gradient(low = "#FFFFFF", high = "#984EA3", 
                          limits = c(min_exp_sen, max_exp_sen)) +
      scale_y_continuous(limits = c(0, 2.9)) +
      xlab("") +
      ylab("") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.45, "cm"))
    if (save_plots) {
      ggsave(paste("figs/fig4/fig4_ts_beninca_2009_",
                   training_prop, "_training_", steps_ahead, "_steps_ahead",
                   "_picophytoplankton.pdf", sep = ""), 
             fig_picophytoplankton, width = 34, height = 5, units = "cm")
    }
    # time series plot for rotifers
    fig_rotifers <- ggplot(ts_plot, aes(x = day, y = rotifers, fill = exp_sens_rotifers)) +
      geom_ribbon(aes(xmin = day[1], xmax = day[start_last_train]),
                  color = "gray90", fill = "gray90") +
      geom_line(size = 0.5) + 
      geom_point(size = 3, shape = 21) + 
      scale_fill_gradient(low = "#FFFFFF", high = "#F781BF", 
                          limits = c(min_exp_sen, max_exp_sen)) +
      scale_y_continuous(limits = c(0, 2.9)) +
      xlab("") +
      ylab("") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.45, "cm"))
    if (save_plots) {
      ggsave(paste("figs/fig4/fig4_ts_beninca_2009_",
                   training_prop, "_training_", steps_ahead, "_steps_ahead",
                   "_rotifers.pdf", sep = ""), 
             fig_rotifers, width = 34, height = 5, units = "cm")
    }
  }
  if (data[i] == "beninca_2015") {
    # adding dates
    beninca_2015_date <- read.csv("data/empirical_time_series/beninca_2015_date.csv")
    beninca_2015_date <- beninca_2015_date$Date
    beninca_2015_date <- beninca_2015_date[1:(nrow(ts) - steps_ahead)]
    ts_plot$date <- as.yearmon(as.Date(beninca_2015_date, "%d/%m/%y"))
    # time series plot for algae
    fig_algae <- ggplot(ts_plot, aes(x = date, y = algae, fill = exp_sens_algae)) +
      geom_ribbon(aes(xmin = date[1], xmax = date[start_last_train]),
                  color = "gray90", fill = "gray90") +
      geom_line(size = 0.5) + 
      geom_point(size = 3, shape = 21) + 
      scale_x_yearmon(n = 10) + 
      scale_fill_gradient(low = "#FFFFFF", high = "#4DAF4A", 
                          limits = c(min_exp_sen, max_exp_sen)) +
      scale_y_continuous(limits = c(0, 95.6)) +
      xlab("") +
      ylab("") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.45, "cm"))
    if (save_plots) {
      ggsave(paste("figs/fig4/fig4_ts_beninca_2015_",
                   training_prop, "_training_", steps_ahead, "_steps_ahead",
                   "_algae.pdf", sep = ""), 
             fig_algae, width = 34, height = 5, units = "cm")
    }
    # time series plot for barnacles
    fig_barnacles <- ggplot(ts_plot, aes(x = date, y = barnacles, fill = exp_sens_barnacles)) +
      geom_ribbon(aes(xmin = date[1], xmax = date[start_last_train]),
                  color = "gray90", fill = "gray90") +
      geom_line(size = 0.5) + 
      geom_point(size = 3, shape = 21) + 
      scale_x_yearmon(n = 10) + 
      scale_fill_gradient(low = "#FFFFFF", high = "#377EB8", 
                          limits = c(min_exp_sen, max_exp_sen)) +
      scale_y_continuous(limits = c(0, 95.6)) +
      xlab("") +
      ylab("") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.45, "cm"))
    if (save_plots) {
      ggsave(paste("figs/fig4/fig4_ts_beninca_2015_",
                   training_prop, "_training_", steps_ahead, "_steps_ahead",
                   "_barnacles.pdf", sep = ""), 
             fig_barnacles, width = 34, height = 5, units = "cm")
    }
    # time series plot for mussels
    fig_mussels <- ggplot(ts_plot, aes(x = date, y = mussels, fill = exp_sens_mussels)) +
      geom_ribbon(aes(xmin = date[1], xmax = date[start_last_train]),
                  color = "gray90", fill = "gray90") +
      geom_line(size = 0.5) + 
      geom_point(size = 3, shape = 21) + 
      scale_x_yearmon(n = 10) + 
      scale_fill_gradient(low = "#FFFFFF", high = "#984EA3", 
                          limits = c(min_exp_sen, max_exp_sen)) +
      scale_y_continuous(limits = c(0, 95.6)) +
      xlab("") +
      ylab("") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.45, "cm"))
    if (save_plots) {
      ggsave(paste("figs/fig4/fig4_ts_beninca_2015_",
                   training_prop, "_training_", steps_ahead, "_steps_ahead",
                   "_mussels.pdf", sep = ""), 
             fig_mussels, width = 34, height = 5, units = "cm")
    }
    # time series plot for rock
    fig_rock <- ggplot(ts_plot, aes(x = date, y = rock, fill = exp_sens_rock)) +
      geom_ribbon(aes(xmin = date[1], xmax = date[start_last_train]),
                  color = "gray90", fill = "gray90") +
      geom_line(size = 0.5) + 
      geom_point(size = 3, shape = 21) + 
      scale_x_yearmon(n = 10) + 
      scale_fill_gradient(low = "#FFFFFF", high = "#F781BF", 
                          limits = c(min_exp_sen, max_exp_sen)) +
      scale_y_continuous(limits = c(0, 95.6)) +
      xlab("") +
      ylab("") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.45, "cm"))
    if (save_plots) {
      ggsave(paste("figs/fig4/fig4_ts_beninca_2015_",
                   training_prop, "_training_", steps_ahead, "_steps_ahead",
                   "_rock.pdf", sep = ""), 
             fig_rock, width = 34, height = 5, units = "cm")
    }
  }
}
