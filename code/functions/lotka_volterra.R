# Returns the rates of change for the Lotka-Volterra model (to use with ode function)

# Arguments:
# t: vector of time steps
# state: vector of initial abundance values
# parameters: vector of parameters (r_i and a_ij values)

# Output:
# dx_dt: abundance rates of change 

lotka_volterra <- function(t, state, parameters) {
  r <- parameters[1:length(state)]
  A <- matrix(parameters[(length(state) + 1):(length(state) + length(state) * length(state))],
              nrow = length(state), ncol = length(state))
  # rate of change
  dx_dt <- state * (r + c(A %*% state))
  # return the rate of change
  list(dx_dt)
}

