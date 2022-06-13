# Code for Fig S7: comparison between expected sensitivity, eigenvector, rate of change, 
# and abundance ranking approaches for multiple models (computes rank correlations by 
# subsetting data using the leading eigenvalue of the Jacobian matrix)

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
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}

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
# whether to save plots
save_plots <- TRUE

# compute correlation between species sensitivities and different rankings ------------------------------
# objects to store results
df_predictions_full <- data.frame()
# loop over models
for (i in 1:length(func_name)) {
  # load results
  load(file = paste("results/perturbation_analyses/jacobian_sensitivities_", func_name[i], "_", n_sp[i], "_sp_",
                    pert_dist, "_pert_", pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, "_",
                    jacobian, "_J_", noise, "_noise", ".RData", sep = ""))
  
  # compute correlation for each eigenavalue percentile ------------------------------ 
  # percentiles of leading eigenvalue to use
  lambda1_perc <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  # data frame to store results
  df_predictions <- data.frame()
  # loop over percentiles
  for (j in 1:length(lambda1_perc)) {
    # obtain subset of data frame
    sub_df <- df[df$lambda_1 > quantile(df$lambda_1, probs = lambda1_perc[j]), ]
    # model name
    model <- data.frame(model = rep(paste(func_name[i], n_sp[i], "sp", sep = "_"), nrow(sub_df)))
    # data frames with observed sensitivity and prediction of each ranking approach
    df_expect_sens <- cbind(sub_df[ , paste("s", 1:n_sp[i], sep = "")], 
                            sub_df[ , paste("expected_sensitivity_", 1:n_sp[i], sep = "")])
    df_eigen_align <- cbind(sub_df[ , paste("s", 1:n_sp[i], sep = "")], 
                            sub_df[ , paste("eigen_alignment_", 1:n_sp[i], sep = "")])
    df_rate_change <- cbind(sub_df[ , paste("s", 1:n_sp[i], sep = "")], 
                            sub_df[ , paste("delta_x", 1:n_sp[i], sep = "")])
    df_abundance <- cbind(sub_df[ , paste("s", 1:n_sp[i], sep = "")], 
                          -sub_df[ , paste("x", 1:n_sp[i], sep = "")])
    # compute rank correlations for each point in time
    predictions_expect_sens <- apply(df_expect_sens, 1, 
                                     function(x) cor(as.numeric(x[1:n_sp[i]]), 
                                                     as.numeric(x[(n_sp[i]+1):(2*n_sp[i])]), 
                                                     method = "spearman"))
    predictions_eigen_align <- apply(df_eigen_align, 1, 
                                     function(x) cor(as.numeric(x[1:n_sp[i]]), 
                                                     as.numeric(x[(n_sp[i]+1):(2*n_sp[i])]), 
                                                     method = "spearman"))
    predictions_rate_change <- apply(df_rate_change, 1, 
                                     function(x) cor(as.numeric(x[1:n_sp[i]]), 
                                                     as.numeric(x[(n_sp[i]+1):(2*n_sp[i])]), 
                                                     method = "spearman"))
    predictions_abundance <- apply(df_abundance, 1, 
                                   function(x) cor(as.numeric(x[1:n_sp[i]]), 
                                                   as.numeric(x[(n_sp[i]+1):(2*n_sp[i])]), 
                                                   method = "spearman"))
    # merge all results into a single data frame 
    df_curr <- cbind(model, predictions_expect_sens, predictions_eigen_align, 
                     predictions_rate_change, predictions_abundance)
    # add eigenvalue percentile information
    df_curr$lambda1_perc <- lambda1_perc[j]
    # merge results from current model to full data frame
    df_predictions <- rbind(df_predictions, df_curr)
  }
  # merge results into single data frame
  df_predictions_full <- rbind(df_predictions_full, df_predictions)
}

# figure showing correlations for each ranking method and each eigenvalue percentile ------------------------------
# computing mean rank correlation per model and eigenvalue percentile
summ_df <- ddply(df_predictions_full, c("model", "lambda1_perc"), summarise,
                 mean_cor_expect_sens = mean(predictions_expect_sens, na.rm = TRUE),
                 mean_cor_eigen_align = mean(predictions_eigen_align, na.rm = TRUE),
                 mean_cor_rate_change = mean(predictions_rate_change, na.rm = TRUE),
                 mean_cor_abundance = mean(predictions_abundance, na.rm = TRUE))
# modifying data frame for plotting
plot_df <- gather(summ_df, "method", "mean_rank_corr", mean_cor_expect_sens:mean_cor_abundance)
# changing model and method names
plot_df$model[plot_df$model == "predator_prey_2_sp"] <- "Predator-prey (2 sp)"
plot_df$model[plot_df$model == "food_chain_3_sp"] <- "Food chain (3 sp)"
plot_df$model[plot_df$model == "lotka_volterra_3_sp"] <- "Food web (3 sp)"
plot_df$model[plot_df$model == "lotka_volterra_4_sp"] <- "Competitors (4 sp)"
plot_df$model[plot_df$model == "consumer_resource_5_sp"] <- "Food web (5 sp)"
plot_df$model <- factor(plot_df$model, levels = c("Predator-prey (2 sp)", "Food chain (3 sp)", "Food web (3 sp)",
                                                  "Competitors (4 sp)", "Food web (5 sp)"))
plot_df$method[plot_df$method == "mean_cor_expect_sens"] <- "Expected sensitivity"
plot_df$method[plot_df$method == "mean_cor_eigen_align"] <- "Eigenvector"
plot_df$method[plot_df$method == "mean_cor_rate_change"] <- "Rate of change"
plot_df$method[plot_df$method == "mean_cor_abundance"] <- "Abundance"
plot_df$method <- factor(plot_df$method, levels = c("Expected sensitivity", "Eigenvector", "Rate of change", "Abundance"))
plot_df$lambda1_perc <- as.factor(plot_df$lambda1_perc)
# plot of mean correlation per model as a function of lambda percentile
fig <- ggplot() +
  geom_point(data = plot_df, aes(x = lambda1_perc, y = mean_rank_corr, 
                                 fill = method, shape = method), size = 3.5) +
  scale_fill_simpsons(name = "Ranking method", 
                      labels = c("Expected sensitivity", "Eigenvector", "Rate of change", "Abundance")) +
  scale_shape_manual(values = c(21, 24, 22, 23),
                     name = "Ranking method", 
                     labels = c("Expected sensitivity", "Eigenvector", "Rate of change", "Abundance")) +
  geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed") +
  xlab(TeX("$\\lambda_1$ percentile")) +
  ylab(expression(atop("Mean correlation" ~ (bar(rho)) ~ "between species",
                       "rankings and sensitivities"))) +
  facet_wrap(~model, nrow = length(func_name), scales = "free_y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 21),
        axis.text.x = element_text(size = 16),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 16),
        legend.key.height = unit(0.8, "cm"),
        legend.key.width = unit(0.8, "cm"))
if (save_plots) {
  ggsave(paste("figs/fig_mean_rank_cor_lambda1_perc_", jacobian, "_",
               pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
               k_type, "_time_step_k_", k, ".pdf", sep = ""), 
         fig, width = 22, height = 30, units = "cm")
}
