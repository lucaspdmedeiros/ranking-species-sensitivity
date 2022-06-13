# Returns the rates of change for the chaotic food chain model in
# Upadhyay 2000 Mathematical and Computer Modelling (to use with ode function)

# Arguments:
# t: vector of time steps
# state: vector of initial abundance values
# parameters: vector of parameters 

# Output:
# dx1: producer rate of change 
# dx2: primary consumer rate of change 
# dx3: secondary consumer rate of change 

food_chain <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    dx1 <- r * x1 * (1 - (x1 / k)) - x1 * ((a1 * x2) / (1 + b1 * x1))
    dx2 <- -s * x2 + h * x1 * x2 - x2 * ((a2 * x3) / (1 + b2 * x2))
    dx3 <- -l * x3 + n * x2 * x3
    # return the rate of change
    list(c(dx1, dx2, dx3))
  })
}