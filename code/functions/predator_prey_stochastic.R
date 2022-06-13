# Solves the stochastic predator-prey model in Yodzis 1989 
# Introduction to theoretical ecology

# Arguments:
# N: vector of initial abundances
# p: vector of parameters
# s: vector of noise magnitudes
# times: vector of integration times

# Output:
# out: data frame with species abundances over time

predator_prey_stochastic <- function(N, p, s, times) {
  # set up diffeqr
  de <- diffeqr::diffeq_setup()
  # define deterministic part
  f <- function(u, p, t) {
    du1 <- p[1] * u[1] * ((p[2] * u[2]^2) / (1 + p[2] * p[3] * u[2]^2)) - p[4] * u[1]
    du2 <- p[5] * u[2] * (1 - (u[2] / p[6])) - u[1] * ((p[2] * u[2]^2) / (1 + p[2] * p[3] * u[2]^2))
    deterministic <- c(du1, du2)
    return(deterministic)
  }
  # define stochastic part
  g <- function(u, p, t) {
    du1 <- p[7] * u[1]
    du2 <- p[8] * u[2]
    stochastic <- c(du1, du2)
    return(stochastic)
  }
  # parameters
  p <- as.numeric(c(p, s))
  # integration time steps
  saveat <- times[2] - times[1]
  # time range
  tspan = as.numeric(range(times))
  # initial condition
  u0 <- as.numeric(N)
  # numerical integration
  prob <- de$SDEProblem(f, g, u0, tspan, p)
  sol <- de$solve(prob, saveat = saveat)
  # organize and return data frame
  out <- as.data.frame(cbind(sol$t, as.data.frame(t(sapply(sol$u, identity)))))
  names(out) <- c("time", paste("x", 1:length(N), sep = ""))
  return(out)
}
