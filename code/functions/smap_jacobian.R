# Function that performs the S-map to infer the Jacobian matrix at different
# points in time using a multivariate time series and returns the Jacobian matrices
# as well as the RMSE for leave one out cross validation

# Arguments:
# ts: time series with time in first column and species abundances in other columns
# theta: kernel parameter

# Output list elements:
# smap_J: Jacobian matrices for each point in time
# rmse: average root mean squared error for each species
# smap_intercept: intercepts of the regression fit

smap_jacobian <- function(ts, theta) {
  # number of species
  n_sp <- ncol(ts) - 1  
  # data points to use
  lib <- nrow(ts)
  # species to use
  cols <- paste("x", 1:n_sp, sep = "")
  # perform s-map for each target species (effect of all species on target species)
  smap_J <- lapply(rep(NA, lib), matrix, nrow = n_sp, ncol = n_sp)
  smap_intercept <- matrix(rep(NA, n_sp * lib), nrow = lib, ncol = n_sp)
  rmse <- rep(NA, n_sp)
  for (i in 1:length(cols)) {
    # target species 
    target <- cols[i]
    # performing s-map using the function block_lnlp from the rEDM package
    smap_output <- block_lnlp(block = ts, method = "s-map", 
                              columns = cols,  target_column = target, theta = theta, 
                              stats_only = FALSE, first_column_time = TRUE, 
                              save_smap_coefficients = TRUE, silent = TRUE)
    # RMSE for this target variable
    rmse[i] <- smap_output$stats$rmse[[1]]
    # s-map coefficients (fitted local regression coefficients through time)
    smap_coeffs <- smap_output$smap_coefficients[[1]]
    # s-map intercept values
    smap_intercept[ , i] <- smap_output$smap_coefficients[[1]]$C0[-1]
    # fill time-varying jacobians 
    for (j in 1:lib) {
      smap_J[[j]][i, ] <- as.numeric(smap_coeffs[j+1, 3:(n_sp+2)])
    }
  }
  return(list(smap_J, rmse, smap_intercept))
}
