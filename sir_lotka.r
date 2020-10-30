require(rstan)
library(tidyverse)
require(bayesplot)
require(ggplot2)
require(loo)
require(stats)
require(posterior)
set.seed(123)

lynx_hare_df <-
  read.csv("hudson-bay-lynx-hare.csv",
           comment.char="#")

model = stan_model("stan/sir_lotka.stan")
expose_stan_functions(model)

N <- length(lynx_hare_df$Year) - 1
ts <- 1:N
y_init <- c(lynx_hare_df$Hare[1], lynx_hare_df$Lynx[1])
y <- as.matrix(lynx_hare_df[2:(N + 1), 2:3])
y <- cbind(y[ , 2], y[ , 1]); # hare, lynx order
lynx_hare_data <- list(N = N, ts = ts, y_init = y_init, y = y)

initial_conditions = c(5, 6)
theta_true = c(0.55, 0.03, 0.8, 0.24)

# Solve the SIR system and format result array
# - theta = [beta, gamma]
solve_sir = function(theta, dt) {
  stan_solve_sir(initial_conditions, ts, theta[1], theta[2], theta[3], theta[4], dt)
}

plot_sir = function(y) {
  y %>%
    as_tibble() %>%
    setNames(c("rabbits", "wolves")) %>%
    mutate(Day = ts) %>%
    gather(which, count, -Day) %>%
    ggplot(aes(Day, count)) +
    geom_line(aes(group = which, color = which)) +
    geom_point(aes(color = which))
}

plot_sir(solve_sir(theta_true, 0.5))

check_reliability = function(theta, dt) {
  y_hat = solve_sir(theta, dt)
  y_hat2 = solve_sir(theta, dt / 2)
  max_abs_err = max(abs(y_hat - y_hat2))
  return(max_abs_err)
}

errors = c()
dts = exp(seq(log(2), log(1e-4), length = 10))
for (dt in dts) {
  errors = c(errors, check_reliability(theta_true, dt))
}

qplot(dts, errors, geom = c("point", "line")) +
  scale_x_log10() +
  scale_y_log10()

dt = 1e-2

dispersion = 5 # noise parameter for negative binomial
mu = solve_sir(theta_true, dt)
y = cbind(rlnorm(length(ts), meanlog = log(mu[, 1]), sdlog = .25),
          rlnorm(length(ts), meanlog = log(mu[, 2]), sdlog = .25))

bind_rows(tibble(t = ts, mu = mu[, 1], y = y[, 1], which = "rabbits"),
          tibble(t = ts, mu = mu[, 2], y = y[, 2], which = "wolves")) %>%
  ggplot() +
  geom_line(aes(t, mu, color = which, group = which)) +
  geom_point(aes(t, y, color = which))

dt_low = 0.5

fit = rstan::sampling(model, list(N = length(ts),
                                  ts = ts,
                                  y = y,
                                  initial_conditions = initial_conditions,
                                  dt = dt_low), cores = 4)

mcmc_trace(fit, c("alpha", "beta", "gamma", "delta"))
