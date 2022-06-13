# Solves the ODE composed of original dynamics and tangent dynamics (variational 
# equation) and returns the full set of state variables over time

# Arguments:
# t: vector of time steps
# state: vector of initial abundance values
# parms: vector of parameters 
# original_func: function returning rate of change of original dynamics
# n_sp: number of state variables

# Output list elements:
# dx_dt: rate of change of original state variables
# dp_dt: rate of change of perturbations (linear tangent dynamics)

original_and_tangent <- function(t, state, parms, original_func, n_sp) {
  # state variables for original dynamics
  state_original <- state[1:n_sp]
  # rate of change of original ODE
  dx_dt <- unlist(original_func(state = state_original, parameters = parms))
  names(dx_dt) <- names(state_original)
  # state variables for tangent dynamics
  state_tangent <- state[(n_sp+1):(2*n_sp)]
  # rate of change of tangent dynamics 
  J <- jacobian.full(y = state_original, func = func, parms = parms)
  dp_dt <- c(J %*% state_tangent)
  names(dp_dt) <- paste("delta", 1:n_sp, sep = "")
  # return the rate of change of all state variables
  list(c(dx_dt, dp_dt))
}
