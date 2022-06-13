# Solves the stochastic Lotka-Volterra model

# Arguments:
# N: vector of initial abundances
# p: vector of parameters (r_i and a_ij values)
# s: vector of noise magnitudes
# times: vector of integration times

# Output:
# out: data frame with species abundances over time

lotka_volterra_stochastic <- function(N, p, s, times) {
  # set up diffeqr
  de <- diffeqr::diffeq_setup()
  # define deterministic part
  f <- function(u, p, t) {
    deterministic <- u * (p$r + c(p$A %*% u))
    return(deterministic)
  }
  # define stochastic part
  g <- function(u, p, t) {
    stochastic <- p$s * u
    return(stochastic)
  }
  # parameters
  p <- list(r = r, A = A, s = s)
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
