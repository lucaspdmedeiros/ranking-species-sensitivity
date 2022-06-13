# Code for Figs 4B and C (also Figs S22-S25): correlation between expected sensitivity, 
# eigenvector, rate of change, and abundance approaches and species forecast errors 
# for each empirical time series

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

# computing correlation between different rankings and forecast errors for each time series ------------------------------ 
# data frame to store results
full_results_df <- data.frame()
for (i in 1:length(data)) {
  # load empirical time series
  ts <- read.csv(paste("data/empirical_time_series/", data[i], ".csv", sep = ""), header = TRUE)
  sp_names <- names(ts)[-1]
  n_sp <- ncol(ts) - 1
  # load s-map and forecast results
  if (normalize) {
    load(paste("results/forecast_analyses/sequential_smap_", data[i], "_normalized_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
    load(paste("results/forecast_analyses/sequential_smap_jacobians_", data[i], "_normalized_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
    load(paste("results/forecast_analyses/sequential_forecast_", data[i], "_normalized_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
  } else {
    load(paste("results/forecast_analyses/sequential_smap_", data[i], "_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
    load(paste("results/forecast_analyses/sequential_smap_jacobians_", data[i], "_",
               training_prop, "_training_", steps_ahead, "_steps_ahead",
               ".RData", sep = ""))
    load(paste("results/forecast_analyses/sequential_forecast_", data[i], "_",
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
  # last training point to start with
  start_last_train <- which(ts$time == smap_df$time[1])
  # first test set point
  start_first_test <- which(ts$time == smap_df$time[2])
  
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
  
  # compute correlations for different percentiles of leading eigenvalue  ------------------------------ 
  # percentiles of leading eigenvalue to use
  lambda1_perc <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  # data frame to store results
  results_df <- data.frame()
  # loop over percentiles
  for (j in 1:length(lambda1_perc)) {
    # obtain subset of data frames
    curr_ts <- ts[start_last_train:(nrow(ts) - steps_ahead), ]
    curr_ts <- curr_ts[smap_df$lambda1 > quantile(smap_df$lambda1, probs = lambda1_perc[j]), ]
    curr_rmse_ratio_df <- rmse_ratio_df[smap_df$lambda1 > quantile(smap_df$lambda1, probs = lambda1_perc[j]), ]
    curr_smap_df <- smap_df[smap_df$lambda1 > quantile(smap_df$lambda1, probs = lambda1_perc[j]), ]
    curr_exp_sens_df <- exp_sens_df[smap_df$lambda1 > quantile(smap_df$lambda1, probs = lambda1_perc[j]), ]
    # vectors to store correlations
    cor_rmse_exp_sens <- c()
    cor_rmse_eigenvector <- c()
    cor_rmse_rate_change <- c()
    cor_rmse_abundance <- c()
    # compute correlations
    for (m in 1:nrow(curr_rmse_ratio_df)) {
      # correlation between RMSE and species expected sensitivities
      cor_rmse_exp_sens[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                                  as.numeric(curr_exp_sens_df[m, ]),
                                  method = "spearman")
      # correlation between RMSE and species eigenvector alignment
      cor_rmse_eigenvector[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                                     abs(as.numeric(curr_smap_df[m, 2:(n_sp+1)])),
                                     method = "spearman")
      # correlation between RMSE and species rate of change
      rate_change <- as.numeric(abs((curr_ts[m + 1, 2:(n_sp+1)] - curr_ts[m, 2:(n_sp+1)]) / 
                                      curr_ts[m, 2:(n_sp+1)]))
      rate_change[is.na(rate_change)] <- 0
      cor_rmse_rate_change[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), rate_change,
                                     method = "spearman")
      # correlation between RMSE and species abundance
      cor_rmse_abundance[m] <- cor(as.numeric(curr_rmse_ratio_df[m, ]), 
                                   -as.numeric(curr_ts[m, 2:(n_sp+1)]),
                                   method = "spearman")
    }
    # data frame for current lambda 1 percentile
    curr_df <- data.frame(cor_rmse_ranking = c(cor_rmse_exp_sens,
                                               cor_rmse_eigenvector,
                                               cor_rmse_rate_change,
                                               cor_rmse_abundance),
                          method = c(rep("Expected sensitivity", length(cor_rmse_exp_sens)),
                                     rep("Eigenvector", length(cor_rmse_eigenvector)),
                                     rep("Rate of change", length(cor_rmse_rate_change)),
                                     rep("Abundance", length(cor_rmse_abundance))),
                          lambda_perc = lambda1_perc[j],
                          data = data[i])
    results_df <- rbind(results_df, curr_df)
  }
  # merge data frames
  full_results_df <- rbind(full_results_df, results_df)  
}

# plot of correlations of different ranking methods for full data ------------------------------ 
# extracting and modifying data
results_df_lambda_perc_0 <- subset(full_results_df, lambda_perc == 0)
results_df_lambda_perc_0$data[results_df_lambda_perc_0$data == "beninca_2009"] <- "Marine plankton community"
results_df_lambda_perc_0$data[results_df_lambda_perc_0$data == "beninca_2015"] <- "Rocky intertidal community"
results_df_lambda_perc_0$data <- factor(results_df_lambda_perc_0$data, levels = c("Rocky intertidal community", "Marine plankton community"))
results_df_lambda_perc_0$method <- factor(results_df_lambda_perc_0$method, levels = c("Expected sensitivity", "Eigenvector", "Rate of change", "Abundance"))
# results of different ranking methods
predictions_exp_sens <- results_df_lambda_perc_0[results_df_lambda_perc_0$method == "Expected sensitivity", ]
predictions_eigen_align <- results_df_lambda_perc_0[results_df_lambda_perc_0$method == "Eigenvector", ]
predictions_rate_change <- results_df_lambda_perc_0[results_df_lambda_perc_0$method == "Rate of change", ]
predictions_abundance <- results_df_lambda_perc_0[results_df_lambda_perc_0$method == "Abundance", ]
# compute fraction of points for which each correlation value is observed (expected sensitivity method)
predictions_exp_sens_table <- lapply(tapply(predictions_exp_sens$cor_rmse_ranking, 
                                                     predictions_exp_sens$data, table), function(x) x / sum(x))
predictions_exp_sens_frac <- as.numeric(unlist(predictions_exp_sens_table))
predictions_exp_sens_corr <- as.numeric(unlist(sapply(predictions_exp_sens_table, names)))
predictions_exp_sens_n <- as.numeric(sapply(predictions_exp_sens_table, length))
# compute fraction of points for which each correlation value is observed (eigenvector method)
predictions_eigen_align_table <- lapply(tapply(predictions_eigen_align$cor_rmse_ranking, 
                                                  predictions_eigen_align$data, table), function(x) x / sum(x))
predictions_eigen_align_frac <- as.numeric(unlist(predictions_eigen_align_table))
predictions_eigen_align_corr <- as.numeric(unlist(sapply(predictions_eigen_align_table, names)))
predictions_eigen_align_n <- as.numeric(sapply(predictions_eigen_align_table, length))
# compute fraction of points for which each correlation value is observed (rate of change method)
predictions_rate_change_table <- lapply(tapply(predictions_rate_change$cor_rmse_ranking, 
                                                 predictions_rate_change$data, table), function(x) x / sum(x))
predictions_rate_change_frac <- as.numeric(unlist(predictions_rate_change_table))
predictions_rate_change_corr <- as.numeric(unlist(sapply(predictions_rate_change_table, names)))
predictions_rate_change_n <- as.numeric(sapply(predictions_rate_change_table, length))
# compute fraction of points for which each correlation value is observed (abundance method)
predictions_abundance_table <- lapply(tapply(predictions_abundance$cor_rmse_ranking, 
                                                  predictions_abundance$data, table), function(x) x / sum(x))
predictions_abundance_frac <- as.numeric(unlist(predictions_abundance_table))
predictions_abundance_corr <- as.numeric(unlist(sapply(predictions_abundance_table, names)))
predictions_abundance_n <- as.numeric(sapply(predictions_abundance_table, length))
# vector with data sets
datasets <- c("Rocky intertidal community", "Marine plankton community")
# plot data frame
plot_df <- data.frame(method = c(rep("Expected sensitivity", length(predictions_exp_sens_corr)),
                                 rep("Eigenvector", length(predictions_eigen_align_corr)),
                                 rep("Rate of change", length(predictions_rate_change_corr)),
                                 rep("Abundance", length(predictions_abundance_corr))),
                      data = c(rep(datasets, times = predictions_exp_sens_n),
                               rep(datasets, times = predictions_eigen_align_n),
                               rep(datasets, times = predictions_rate_change_n),
                               rep(datasets, times = predictions_abundance_n)),
                      rank_corr = c(predictions_exp_sens_corr, predictions_eigen_align_corr,
                                    predictions_rate_change_corr, predictions_abundance_corr),
                      frac = c(predictions_exp_sens_frac,
                               predictions_eigen_align_frac,
                               predictions_rate_change_frac,
                               predictions_abundance_frac))
plot_df$data <- factor(plot_df$data, levels = c("Rocky intertidal community", "Marine plankton community"))
plot_df$method <- factor(plot_df$method, levels = c("Expected sensitivity", "Eigenvector", "Rate of change", "Abundance"))
# computing mean rank correlation for each method and model
plot_df$sum_corr_frac <- plot_df$rank_corr * plot_df$frac
summ_df <- ddply(plot_df, c("data", "method"), summarise,
                 mean_rank_corr = sum(sum_corr_frac))
# plot
fig <- ggplot() +
  geom_point(data = plot_df, mapping = aes(x = method, y = rank_corr, size = 100*frac),
             color = "gray80") +
  geom_point(data = summ_df, mapping = aes(x = method, y = mean_rank_corr),
             shape = 95, size = 23, show.legend = FALSE) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  facet_wrap(~data, ncol = 5) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(size = guide_legend(title = "Percentage of\npoints (%)")) + 
  xlab("Ranking method") +
  ylab(expression(atop("Correlation" ~ (rho) ~ "between species",
                       "rankings and forecast errors"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 17),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.7, "cm"))
if (save_plots) {
  ggsave(paste("figs/fig4/fig4_cor_ranking_methods_", training_prop, "_training_", steps_ahead, 
               "_steps_ahead_forecast", ".pdf", sep = ""), 
         fig, width = 28, height = 18, units = "cm")
}

# plot of correlations of different ranking methods for lambda_1 percentile subsets of data ------------------------------ 
# plot data frame
plot_df <- ddply(full_results_df, c("data", "lambda_perc", "method"), summarise,
                 mean_cor_rmse_ranking = mean(cor_rmse_ranking, na.rm = TRUE))
plot_df$data[plot_df$data == "beninca_2009"] <- "Marine plankton community"
plot_df$data[plot_df$data == "beninca_2015"] <- "Rocky intertidal community"
plot_df$data <- factor(plot_df$data, levels = c("Rocky intertidal community", "Marine plankton community"))
plot_df$method <- factor(plot_df$method, levels = c("Expected sensitivity", "Eigenvector", "Rate of change", "Abundance"))
plot_df$lambda_perc <- as.factor(100*plot_df$lambda_perc)
# plot
fig <- ggplot() +
  geom_point(data = plot_df, aes(x = method, y = mean_cor_rmse_ranking, color = lambda_perc), 
             shape = 95, size = 23) +
  geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed") +
  xlab("Ranking method") +
  ylab(expression(atop("Mean correlation" ~ (bar(rho)) ~ "between species",
                       "rankings and forecast errors"))) +
  facet_wrap(~data) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_viridis(name = TeX("$\\lambda_1$ percentile"), option = "plasma", discrete = TRUE,
                      guide = guide_legend(direction = "horizontal", title.hjust = 0.03, title.vjust = 0,
                                           title.position = "top", override.aes = list(size = 3.4, shape = 15))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 17),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 14),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.65, "cm"),
        legend.position = c(0.24, 0.12))
if (save_plots) {
  ggsave(paste("figs/fig4/fig4_cor_ranking_methods_lambda1_", training_prop, "_training_", steps_ahead, 
               "_steps_ahead_forecast", ".pdf", sep = ""), 
         fig, width = 24, height = 18, units = "cm")
}
