# Returns the rates of change for the predator-prey with type II 
# functional response model in Yodzis 1989 Introduction to 
# theoretical ecology (to use with ode function)

# Arguments:
# t: vector of time steps
# state: vector of initial abundance values
# parameters: vector of parameters 

# Output:
# dx1: predator rate of change 
# dx2: prey rate of change 

predator_prey <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    dx1 <- k * x1 * ((a * x2^2) / (1 + a * h * x2^2)) - d * x1
    dx2 <- r * x2 * (1 - (x2 / K)) - x1 * ((a * x2^2) / (1 + a * h * x2^2))
    # return the rate of change
    list(c(dx1, dx2))
  })
}