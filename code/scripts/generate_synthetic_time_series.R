# Generate synthetic time series using different models and also generate plots
# for each attractor in state space (Fig S5)

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
source("code/functions/lotka_volterra_stochastic.R")
source("code/functions/predator_prey.R")
source("code/functions/predator_prey_stochastic.R")
source("code/functions/food_chain.R")
source("code/functions/food_chain_stochastic.R")
source("code/functions/consumer_resource.R")
source("code/functions/consumer_resource_stochastic.R")
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(plotly)) {install.packages("plotly"); library(plotly)}
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if(!require(diffeqr)) {install.packages("diffeqr"); library(diffeqr)}
if(!require(JuliaCall)) {install.packages("JuliaCall"); library(JuliaCall)}

# simulate model to generate time series ------------------------------
# set seed to replicate time series
set.seed(42)
# whether to save plots
save_plots <- TRUE
# number of species
n_sp <- 2
# model to use
func_name <- "predator_prey"
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# amount of observational noise to add to time series
noise <- 0
# load model settings
source("code/scripts/model_settings.R")
# perform numerical integration
if (grepl("stochastic", func_name)) {
  ts <- func(N = state, p = parms, s = s, times = times)
} else {
  ts <- as.data.frame(ode(y = state, times = times, func = func, parms = parms, method = "ode45"))
}
# remove transient
ts <- tail(ts, ts_length * sampling_freq)
# sample points according to sampling frequency
ts <- ts[seq(1, nrow(ts), by = sampling_freq), ]
# simplifying times
ts$time <- 1:nrow(ts)

# adding observational noise ------------------------------
if (noise > 0) {
  # time series without time column
  ts_abund <- ts[ , -1]
  # adding noise
  ts[ , 2:(n_sp+1)] <- apply(ts_abund, 2, function(x, noise) {
    x + rnorm(length(x), 0, noise * x)
  }, noise = noise)
  # fixing potential negative abundances
  ts[ts < 0] <- 0.00001
}

# plotting attractor ------------------------------
if (n_sp == 2) {
  # plot for 2 species
  fig <- plot_ly(data = ts, x = ~x1, y = ~x2,
                 type = 'scatter', mode = "lines+markers",
                 line = list(width = 2.5, color = "black"),
                 marker = list(size = 11, color = "black")) %>% 
    layout(xaxis = list(title = "N<sub>1</sub>",
                        titlefont = list(size = 42, 
                                         family = "Arial, sans-serif"), 
                        tickfont = list(size = 30,
                                        family = "Arial, sans-serif"),
                        ticklen = 4,
                        tickcolor = "#7F7F7F",
                        gridcolor = "#7F7F7F",
                        gridwidth = 0.8,
                        zerolinewidth = 0),
           yaxis = list(title = "N<sub>2</sub>",
                        titlefont = list(size = 42, 
                                         family = "Arial, sans-serif"),
                        tickfont = list(size = 30,
                                        family = "Arial, sans-serif"),
                        ticklen = 4,
                        tickcolor = "#7F7F7F",
                        gridcolor = "#7F7F7F",
                        gridwidth = 0.8,
                        zerolinewidth = 0),
           margin = list(l = 0,
                         r = 0,
                         b = 0,
                         t = 0))
  # save plot
  if (save_plots) {
    orca(fig, paste("figs/figSI/figS5/figSI_attractor_", func_name, "_", n_sp, "_sp_",
                    noise, "_noise", ".pdf", sep = ""), format = "pdf",
         width = 600, height = 600)
  }
} else {
  # plot for more than 2 species
  fig <- plot_ly(data = ts, x = ~x1, y = ~x2, z = ~x3, 
                 type = 'scatter3d', mode = "lines+markers",
                 line = list(width = 2, color = "black"),
                 marker = list(size = 4.5, color = "black")) %>% 
    layout(scene = list(xaxis = list(title = "N<sub>1</sub>",
                                     titlefont = list(size = 30, 
                                                      family = "Arial, sans-serif"), 
                                     tickfont = list(size = 16,
                                                     family = "Arial, sans-serif"),
                                     ticklen = 4,
                                     tickcolor = "#7F7F7F",
                                     gridcolor = "#7F7F7F",
                                     gridwidth = 0.6,
                                     zerolinewidth = 0),
                        yaxis = list(title = "N<sub>2</sub>",
                                     titlefont = list(size = 30, 
                                                      family = "Arial, sans-serif"),
                                     tickfont = list(size = 16,
                                                     family = "Arial, sans-serif"),
                                     ticklen = 4,
                                     tickcolor = "#7F7F7F",
                                     gridcolor = "#7F7F7F",
                                     gridwidth = 0.6,
                                     zerolinewidth = 0),
                        zaxis = list(title = "N<sub>3</sub>",
                                     titlefont = list(size = 30, 
                                                      family = "Arial, sans-serif"),
                                     tickfont = list(size = 16,
                                                     family = "Arial, sans-serif"),
                                     ticklen = 4,
                                     tickcolor = "#7F7F7F",
                                     gridcolor = "#7F7F7F",
                                     gridwidth = 0.6,
                                     zerolinewidth = 0),
                        camera = list(eye = list(x = 1.8, y = 0.8, z = 1)),
                        aspectratio = list(x = 0.8, y = 0.8, z = 0.8)),
           margin = list(l = 0,
                         r = 0,
                         b = 0,
                         t = 0))
  # save plot
  if (save_plots) {
    orca(fig, paste("figs/figSI/figS5/figSI_attractor_", func_name, "_", n_sp, "_sp_",
                    noise, "_noise", ".pdf", sep = ""), format = "pdf",
         width = 1000, height = 600)
  }
}

# saving object as RData ------------------------------
save(ts, file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                      ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                      ".RData", sep = ""))
