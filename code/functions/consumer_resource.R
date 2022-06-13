# Returns the rates of change for the chaotic food web model in
# Deyle et al 2016 (to use with ode function)

# Arguments:
# t: vector of time steps
# state: vector of initial abundance values
# parameters: vector of parameters 

# Output:
# dx1: secondary consumer rate of change 
# dx2: secondary consumer rate of change 
# dx3: primary consumer rate of change 
# dx4: primary consumer rate of change 
# dx5: producer rate of change 

consumer_resource <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    dx1 <- v1 * l1 * ((x1 * x3) / (x3 + x3_par)) - v1 * x1
    dx2 <- v2 * l2 * ((x2 * x4) / (x4 + x4_par)) - v2 * x2
    dx3 <- u1 * k1 * ((x3 * x5) / (x5 + x5_par)) - v1 * l1 * ((x1 * x3) / (x3 + x3_par)) - u1 * x3
    dx4 <- u2 * k2 * ((x4 * x5) / (x5 + x5_par)) - v2 * l2 * ((x2 * x4) / (x4 + x4_par)) - u2 * x4
    dx5 <- x5 * (1 - (x5 / k)) - (u1 * k1 * ((x3 * x5) / (x5 + x5_par)) + u2 * k2 * ((x4 * x5) / (x5 + x5_par)))
    # return the rate of change
    list(c(dx1, dx2, dx3, dx4, dx5))
  })
}