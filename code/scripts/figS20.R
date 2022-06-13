# Code for Fig S20: forecast errors, expected sensitivities, and
# eigenvector alignments for each species over time for each
# empirical time series

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
k <- 3

# compute and plot expected sensitivities ------------------------------ 
for (i in 1:length(data)) {
  # load empirical time series
  ts <- read.csv(paste("data/empirical_time_series/", data[i], ".csv", sep = ""), header = TRUE)
  sp_names <- names(ts)[-1]
  n_sp <- ncol(ts) - 1
  # load s-map results
  load(paste("results/forecast_analyses/sequential_smap_", data[i], "_",
             training_prop, "_training_", steps_ahead, "_steps_ahead",
             ".RData", sep = ""))
  load(paste("results/forecast_analyses/sequential_smap_jacobians_", data[i], "_",
             training_prop, "_training_", steps_ahead, "_steps_ahead",
             ".RData", sep = ""))
  # load forecast results
  load(paste("results/forecast_analyses/sequential_forecast_", data[i], "_",
             training_prop, "_training_", steps_ahead, "_steps_ahead",
             ".RData", sep = ""))
  
  # computing expected species sensitivities ------------------------------ 
  # last training point to start with
  start_last_train <- which(ts$time == smap_df$time[1])
  # first test set point
  start_first_test <- which(ts$time == smap_df$time[2])
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
  # change names
  names(exp_sens_df) <- sp_names
  # scaling expected sensitivities
  exp_sens_df <- exp_sens_df / apply(exp_sens_df, 1, sum)
  # adding time to data frame
  exp_sens_df$time <- ts$time[start_last_train:(nrow(ts) - steps_ahead)]
  # adding variable name
  exp_sens_df$type <- "Expected sensitivity"
  
  # computing species eigenvector alignments ------------------------------ 
  # create data frame with eigenvector alignments
  eigen_align_df <- abs(smap_df[ , paste("v1", sp_names, sep = "_")])
  # change names
  names(eigen_align_df) <- sp_names
  # scaling eigenvector alignments
  eigen_align_df <- eigen_align_df / apply(eigen_align_df, 1, sum)
  # adding time to data frame
  eigen_align_df$time <- ts$time[start_last_train:(nrow(ts) - steps_ahead)]
  # adding variable name
  eigen_align_df$type <- "Eigenvector alignment"
  
  # computing forecast errors sequentially ------------------------------ 
  # total number of forecasts
  n_forecasts <- nrow(smap_df)
  # objects to store results
  rmse_ratio <- list()
  # loop to compute RMSE of sequential forecasts
  for (j in 1:n_forecasts) {
    # forecasts for time smap$time[j]
    forecast <- forecast_df[forecast_df$last_train_time == smap_df$time[j], -c(1, 2)]
    # just test set
    ts_test <- ts[(start_first_test + j - 1):(start_first_test + j + steps_ahead - 2), -1]
    # test set standard deviation
    test_sd <- apply(ts_test, 2, sd)
    # add a small value when sd is zero
    test_sd[test_sd == 0] <- 0.001
    # forecast RMSE
    forecast_rmse <- apply((ts_test - forecast)^2, 2, function(x) sqrt(mean(x)))
    # naive forecasts
    naive_forecast <- matrix(rep(as.numeric(ts[ts$time == smap_df$time[j], -1]), each = steps_ahead), 
                             nrow = steps_ahead,
                             ncol = n_sp, byrow = FALSE)
    # add a small value when naive forecast is zero
    naive_forecast[naive_forecast == 0] <- 0.001
    # naive forecast RMSE
    naive_rmse <- apply((ts_test - naive_forecast)^2, 2, function(x) sqrt(mean(x)))
    # standardizing forecast RMSE
    rmse_ratio[[j]] <- forecast_rmse / naive_rmse
  }
  # transforming list to data frame
  rmse_ratio_df <- as.data.frame(do.call("rbind", rmse_ratio))
  # change names
  names(rmse_ratio_df) <- sp_names
  # scaling forecast errors
  rmse_ratio_df <- rmse_ratio_df / apply(rmse_ratio_df, 1, sum)
  # adding time to data frame
  rmse_ratio_df$time <- ts$time[start_last_train:(nrow(ts) - steps_ahead)]
  # adding variable name
  rmse_ratio_df$type <- "Forecast error"
  
  # plotting species sensitivities over time ------------------------------ 
  # merge all data frames 
  full_df <- rbind(exp_sens_df, eigen_align_df, rmse_ratio_df)
  # creating data frame for plotting
  plot_df <- gather(data = full_df, key = "species", value = "variable", -c(time, type))
  plot_df$type <- factor(plot_df$type, levels = c("Forecast error", 
                                                  "Expected sensitivity",
                                                  "Eigenvector alignment"))
  # plot settings
  if (data[i] == "beninca_2009") {
    width <- 3.35
    plot_df$species <- factor(plot_df$species, levels = c("calanoids", "nanophytoplankton",
                                                          "picophytoplankton", "rotifers"))
  }
  if (data[i] == "beninca_2015") {
    width <- 1
    plot_df$species <- factor(plot_df$species, levels = c("rock", "barnacles",
                                                          "algae", "mussels"))
  }
  # plot
  fig <- ggplot(data = plot_df, aes(x = time, y = variable, fill = species)) +
    geom_bar(stat = "identity", width = width) +
    scale_fill_manual(values = palette) +
    facet_wrap(~type, ncol = 3) +
    xlab("Time") +
    ylab("Scaled variable") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          strip.text = element_text(size = 14),
          strip.background = element_rect(fill = "white", size = 1.5),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 11),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.key.size = unit(0.5, "cm"),
          legend.position = "top")
  # save plot
  if (save_plots) {
    ggsave(paste("figs/fig_forecast_errors_rankings_over_time_", data[i], "_",
                 training_prop, "_training_", steps_ahead, "_steps_ahead",
                 ".pdf", sep = ""), 
           fig, width = 30, height = 10, units = "cm")
  }
}
