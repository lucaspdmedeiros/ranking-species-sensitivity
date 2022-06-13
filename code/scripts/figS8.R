# Code for Fig S8: inference quality of leading eigenvalue and eigenvector
# of the Jacobian matrix using the S-map

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
source("code/functions/predator_prey.R")
source("code/functions/food_chain.R")
source("code/functions/consumer_resource.R")
source("code/functions/smap_jacobian.R")
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
if(!require(rEDM)) {install.packages("rEDM"); library(rEDM)}

# models and simulation settings ------------------------------
# number of species
n_sp_vec <- c(2, 3, 3, 4, 5)
# model to use
func_name_vec <- c("predator_prey", "food_chain", "lotka_volterra", 
                   "lotka_volterra", "consumer_resource")
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# amount of observational noise to add to time series (0 or 0.1)
noise <- 0
# whether to normalize data for s-map
normalize <- FALSE
# fraction of time series to train s-map
frac_train <- 0.5
# kernel parameter for s-map
theta_seq <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 
               1, 2, 3, 4, 5, 6, 7, 8)
# whether to save plots
save_plots <- TRUE

# infer Jacobian and its eigenvalues and eigenvectors ------------------------------ 
# data frame to store results
results_df <- data.frame()
# loop over models
for (i in 1:length(func_name_vec)) {
  print(func_name_vec[i])
  # number of species
  n_sp <- n_sp_vec[i]
  # model to use
  func_name <- func_name_vec[i]
  # load model settings
  source("code/scripts/model_settings.R")
  # load time series 
  load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                    ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                    ".RData", sep = ""))
  # compute Jacobian sequentially with the s-map using only past data
  # list to store results
  J_smap <- list()
  # time points to train s-map
  train_points <- 1:ceiling(nrow(ts) * frac_train)
  # time points to infer Jacobians
  test_points <- ceiling(nrow(ts) * frac_train):nrow(ts)
  # time series to train s-map
  ts_train <- ts[train_points, ]
  # time series to infer Jacobians
  ts <- ts[test_points, ]
  # infer Jacobian sequentially
  for (j in 1:nrow(ts)) {
    print(j)
    ts_train_curr <- ts_train
    # normalize time series
    if (normalize) {
      ts_train_curr[ , -1] <- scale(ts_train_curr[ , -1])
    }
    # select kernel parameter using cross validation
    rmse_seq <- sapply(theta_seq, function(x) mean(smap_jacobian(x, ts = ts_train_curr)[[2]]))
    theta <- theta_seq[which.min(rmse_seq)]
    # fit smap to past points
    J_train <- smap_jacobian(ts_train_curr, theta = theta)[[1]]
    # save last jacobian
    J_smap[[j]] <- J_train[[length(J_train)]]
    # add a point to training data
    ts_train <- rbind(ts_train, ts[j + 1, ])
    # remove first point of training data
    ts_train <- ts_train[-1, ]
  }
  # analytical Jacobian matrix for each time step
  J_analytical <- dlply(ts, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                              func = func,
                                                              parms = parms))
  # list with eigendecomposition of J matrices
  eigen_J_analytical <- lapply(J_analytical, eigen)
  eigen_J_smap <- lapply(J_smap, eigen)
  # extract leading eigenvalue and eigenvector
  lambda1_analytical <- c()
  lambda1_smap <- c()
  v1_analytical <- list()
  v1_smap <- list()
  cos_v1 <- c()
  cos_random <- c()
  for (j in 1:length(J_analytical)) {
    # analytical leading eigenvalue and eigenvector
    order_analytical <- order(Re(eigen_J_analytical[[j]]$values), decreasing = TRUE)
    lambda1_analytical[j] <- Re(eigen_J_analytical[[j]]$values)[order_analytical[1]]
    v1_analytical[[j]] <- Re(eigen_J_analytical[[j]]$vectors)[ , order_analytical[1]]
    # s-map leading eigenvalue and eigenvector
    order_smap <- order(Re(eigen_J_smap[[j]]$values), decreasing = TRUE)
    lambda1_smap[j] <- Re(eigen_J_smap[[j]]$values)[order_smap[1]]
    v1_smap[[j]] <- Re(eigen_J_smap[[j]]$vectors)[ , order_smap[1]]
    # alignment between analytical and smap eigenvector
    x <- v1_analytical[[j]]
    y <- v1_smap[[j]]
    x <- x / sqrt(sum(x^2))
    y <- y / sqrt(sum(y^2))
    cos_v1[j] <- abs(sum(x * y))
    # alignment between two random vectors
    x <- rnorm(n_sp, 0, 1)
    x <- x / sqrt(sum(x^2))
    y <- rnorm(n_sp, 0, 1)
    y <- y / sqrt(sum(y^2))
    cos_random[j] <- abs(sum(x * y))
  }
  # build data frame for the results of current model
  curr_df <- data.frame(model = func_name,
                        n_sp = n_sp,
                        time = ts$time,
                        lambda1_analytical = lambda1_analytical,
                        lambda1_smap = lambda1_smap,
                        cos_v1 = cos_v1,
                        cos_random = cos_random)
  # merge with full data frame
  results_df <- rbind(results_df, curr_df)
}

# inference quality of leading eigenvalue ------------------------------ 
# change models in data frame 
results_df$model[results_df$model == "predator_prey"] <- "Predator-prey (2 sp)"
results_df$model[results_df$model == "food_chain"] <- "Food chain (3 sp)"
results_df$model[results_df$model == "lotka_volterra" & results_df$n_sp == 3] <- "Food web (3 sp)"
results_df$model[results_df$model == "lotka_volterra" & results_df$n_sp == 4] <- "Competitors (4 sp)"
results_df$model[results_df$model == "consumer_resource"] <- "Food web (5 sp)"
results_df$model <- factor(results_df$model, levels = c("Predator-prey (2 sp)", "Food chain (3 sp)", "Food web (3 sp)",
                                                        "Competitors (4 sp)", "Food web (5 sp)"))
# compute correlation between true and inferred eigenvalue for each model
ddply(results_df, c("model", "n_sp"), summarise,
      cor_lambda1 = cor(lambda1_analytical, lambda1_smap, method = "pearson"))

# plot of inference quality of leading eigenvector ------------------------------ 
# create new data frame for plotting
plot_df <- gather(results_df, "type", "cos", cos_v1:cos_random)
plot_df$type[plot_df$type == "cos_v1"] <- "True and inferred\neigenvector"
plot_df$type[plot_df$type == "cos_random"] <- "Two random vectors"
plot_df$type <- factor(plot_df$type, levels = c("True and inferred\neigenvector", "Two random vectors"))
# plot
fig <- ggplot(data = plot_df, aes(x = type, y = cos, group = type)) +
  geom_boxplot(size = 1.2, outlier.size = 1.5) +
  facet_wrap(~model, nrow = 3) +
  xlab("") +
  ylab(TeX("Cosine of the angle (alignment) between the two vectors")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white", size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 21),
        legend.position = "none")
if (save_plots) {
  ggsave("figs/fig_smap_inference_quality_eigenvector.pdf", fig, width = 20, height = 30, units = "cm")
}
