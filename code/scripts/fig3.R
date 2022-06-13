# Code for Fig 3 (also Figs S9-S18): comparison between expected sensitivity, 
# eigenvector, rate of change, and abundance ranking approaches for multiple models

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
# whether to save plots
save_plots <- TRUE

# compute correlation between observed sensitivities and different rankings ------------------------------
# objects to store results
df_predictions <- data.frame()
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
  # model name
  model <- data.frame(model = rep(paste(func_name[i], n_sp[i], "sp", sep = "_"), nrow(df)))
  # data frames with observed sensitivity and prediction of each ranking approach
  df_expect_sens <- cbind(df[ , paste("s", 1:n_sp[i], sep = "")], 
                          df[ , paste("expected_sensitivity_", 1:n_sp[i], sep = "")])
  df_eigen_align <- cbind(df[ , paste("s", 1:n_sp[i], sep = "")], 
                          df[ , paste("eigen_alignment_", 1:n_sp[i], sep = "")])
  df_rate_change <- cbind(df[ , paste("s", 1:n_sp[i], sep = "")], 
                          df[ , paste("delta_x", 1:n_sp[i], sep = "")])
  df_abundance <- cbind(df[ , paste("s", 1:n_sp[i], sep = "")], 
                        -df[ , paste("x", 1:n_sp[i], sep = "")])
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
  # merge results from current model to full data frame
  df_predictions <- rbind(df_predictions, df_curr)
}
# change model to factor
df_predictions$model <- factor(df_predictions$model, levels = c("predator_prey_2_sp", "food_chain_3_sp", "lotka_volterra_3_sp",
                                                                "lotka_volterra_4_sp", "consumer_resource_5_sp"))

# figure showing correlations for each ranking method ------------------------------
# compute fraction of points for which each correlation value is observed (expected sensitivity method)
predictions_expect_sens_table <- lapply(tapply(df_predictions$predictions_expect_sens, df_predictions$model, table),
                                        function(x) x / sum(x))
predictions_expect_sens_frac <- as.numeric(unlist(predictions_expect_sens_table))
predictions_expect_sens_corr <- as.numeric(unlist(sapply(predictions_expect_sens_table, names)))
predictions_expect_sens_n <- as.numeric(sapply(predictions_expect_sens_table, length))
# compute fraction of points for which each correlation value is observed (eigenvector method)
predictions_eigen_align_table <- lapply(tapply(df_predictions$predictions_eigen_align, df_predictions$model, table),
                                           function(x) x / sum(x))
predictions_eigen_align_frac <- as.numeric(unlist(predictions_eigen_align_table))
predictions_eigen_align_corr <- as.numeric(unlist(sapply(predictions_eigen_align_table, names)))
predictions_eigen_align_n <- as.numeric(sapply(predictions_eigen_align_table, length))
# compute fraction of points for which each correlation value is observed (rate of change method)
predictions_rate_change_table <- lapply(tapply(df_predictions$predictions_rate_change, df_predictions$model, table),
                                                      function(x) x / sum(x))
predictions_rate_change_frac <- as.numeric(unlist(predictions_rate_change_table))
predictions_rate_change_corr <- as.numeric(unlist(sapply(predictions_rate_change_table, names)))
predictions_rate_change_n <- as.numeric(sapply(predictions_rate_change_table, length))
# compute fraction of points for which each correlation value is observed (abundance method)
predictions_abundance_table <- lapply(tapply(df_predictions$predictions_abundance, df_predictions$model, table),
                                                 function(x) x / sum(x))
predictions_abundance_frac <- as.numeric(unlist(predictions_abundance_table))
predictions_abundance_corr <- as.numeric(unlist(sapply(predictions_abundance_table, names)))
predictions_abundance_n <- as.numeric(sapply(predictions_abundance_table, length))
# change model names
models <- c("Predator-prey (2 sp)", "Food chain (3 sp)", "Food web (3 sp)",
            "Competitors (4 sp)", "Food web (5 sp)")
# create data frame for plotting
plot_df <- data.frame(method = c(rep("Expected sensitivity", length(predictions_expect_sens_corr)),
                                 rep("Eigenvector", length(predictions_eigen_align_corr)),
                                 rep("Rate of change", length(predictions_rate_change_corr)),
                                 rep("Abundance", length(predictions_abundance_corr))),
                      model = c(rep(models, times = predictions_expect_sens_n),
                                rep(models, times = predictions_eigen_align_n),
                                rep(models, times = predictions_rate_change_n),
                                rep(models, times = predictions_abundance_n)),
                      rank_corr = c(predictions_expect_sens_corr, predictions_eigen_align_corr,
                                    predictions_rate_change_corr, predictions_abundance_corr),
                      frac = c(predictions_expect_sens_frac,
                               predictions_eigen_align_frac,
                               predictions_rate_change_frac,
                               predictions_abundance_frac))
plot_df$method <- factor(plot_df$method, levels = c("Expected sensitivity", "Eigenvector", "Rate of change", "Abundance"))
# changing model to factor
plot_df$model <- factor(plot_df$model, levels = c("Predator-prey (2 sp)", "Food chain (3 sp)", "Food web (3 sp)",
                                                  "Competitors (4 sp)", "Food web (5 sp)"))
# computing mean rank correlation for each method and model
plot_df$sum_corr_frac <- plot_df$rank_corr * plot_df$frac
summ_df <- ddply(plot_df, c("model", "method"), summarise,
                 mean_rank_corr = sum(sum_corr_frac))
# plot
fig <- ggplot() +
  geom_point(data = plot_df, mapping = aes(x = method, y = rank_corr, size = 100*frac),
             color = "gray80") +
  geom_point(data = summ_df, mapping = aes(x = method, y = mean_rank_corr),
             shape = 95, size = 22, show.legend = FALSE) +
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  facet_wrap(~model, ncol = 5) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(size = guide_legend(title = "Percentage of\npoints (%)")) + 
  xlab("Ranking method") +
  ylab(expression(atop("Correlation" ~ (rho) ~ "between species",
                       "rankings and sensitivities"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 17),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 17, angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.key.size = unit(0.8, "cm"))
if (save_plots) {
  if (k_error == "none") {
    if (normalize) {
      ggsave(paste("figs/fig3/fig3_", jacobian, "_normalized_", pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                   k_type, "_time_step_k_", k, "_", noise, "_noise", ".pdf", sep = ""), 
             fig, width = 38, height = 16, units = "cm")
    } else {
      ggsave(paste("figs/fig3/fig3_", jacobian, "_", pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                   k_type, "_time_step_k_", k, "_", noise, "_noise", ".pdf", sep = ""), 
             fig, width = 38, height = 16, units = "cm")
    }
  }
  if (k_error == "misspecified") {
    ggsave(paste("figs/fig3/fig3_", jacobian, "_", pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                 k_type, "_time_step_k_", k, "_misspecified_", noise, "_noise", ".pdf", sep = ""), 
           fig, width = 38, height = 16, units = "cm")
  }
  if (k_error == "greatly_misspecified") {
    ggsave(paste("figs/fig3/fig3_", jacobian, "_", pert_dist, "_pert_", pert_magnitude, "_pert_magn_", 
                 k_type, "_time_step_k_", k, "_greatly_misspecified_", noise, "_noise", ".pdf", sep = ""), 
           fig, width = 38, height = 16, units = "cm")
  }
}
