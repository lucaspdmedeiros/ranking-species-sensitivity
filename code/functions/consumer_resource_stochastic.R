# Solves the stochastic food web model in Deyle et al 2016

# Arguments:
# N: vector of initial abundances
# p: vector of parameters
# s: vector of noise magnitudes
# times: vector of integration times

# Output:
# out: data frame with species abundances over time

consumer_resource_stochastic <- function(N, p, s, times) {
  # set up diffeqr
  de <- diffeqr::diffeq_setup()
  # define deterministic part
  f <- function(u, p, t) {
    du1 <- p[1] * p[3] * ((u[1] * u[3]) / (u[3] + p[5])) - p[1] * u[1]
    du2 <- p[2] * p[4] * ((u[2] * u[4]) / (u[4] + p[6])) - p[2] * u[2]
    du3 <- p[7] * p[9] * ((u[3] * u[5]) / (u[5] + p[11])) - p[1] * p[3] * ((u[1] * u[3]) / (u[3] + p[5])) - p[7] * u[3]
    du4 <- p[8] * p[10] * ((u[4] * u[5]) / (u[5] + p[11])) - p[2] * p[4] * ((u[2] * u[4]) / (u[4] + p[6])) - p[8] * u[4]
    du5 <- u[5] * (1 - (u[5] / p[12])) - (p[7] * p[9] * ((u[3] * u[5]) / (u[5] + p[11])) + p[8] * p[10] * ((u[4] * u[5]) / (u[5] + p[11])))
    deterministic <- c(du1, du2, du3, du4, du5)
    return(deterministic)
  }
  # define stochastic part
  g <- function(u, p, t) {
    du1 <- p[13] * u[1]
    du2 <- p[14] * u[2]
    du3 <- p[15] * u[3]
    du4 <- p[16] * u[4]
    du5 <- p[17] * u[5]
    stochastic <- c(du1, du2, du3, du4, du5)
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
