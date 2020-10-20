# Requirements
require(rstan)
library(tidyverse)
require(bayesplot)
require(ggplot2)
require(loo)
require(stats)
require(posterior)
set.seed(123)

# Build model and expose functions
model = stan_model("../stan/sir.stan")
expose_stan_functions(model)


N = 16                                # number of days measured
M = 1000                              # population size
I0 = 20                               # number of infected on day 0
initial_conditions = c(M - I0, I0, 0) # S, I, R on day 0
ts = seq_len(N)                       # measurement times
theta_true = c(1, 0.2)                # true parameter values

# Solve the SIR system and format result array
# - theta = [beta, gamma]
solve_sir = function(theta, rtol, atol, max_num_steps) {
  stan_solve_sir(initial_conditions, ts, theta,
                 c(0.0), M, rtol, atol, max_num_steps) %>%
    unlist %>%
    matrix(ncol = 3, byrow = TRUE)
}

# Function for plotting a solution y_hat
plot_sir = function(y0, y_hat) {
  yyy = rbind(y0, y_hat)
  ttt = c(0.0, ts)
  y = as.vector(yyy)
  day = rep(ttt, 3)
  comp = rep( c("S", "I", "R"), each = length(ttt))
  df = data.frame(day, y, as.factor(comp))
  colnames(df) = c("Day", "y", "Compartment")
  aes = aes_string(x = "Day", y = "y", group = "Compartment",
                   color = "Compartment")
  plt = ggplot(df, aes) + geom_line(lwd = 1) + geom_point()
  return(plt)
}

plot_sir(initial_conditions, solve_sir(theta_true, 1e-4, 1e-4, 1e8))


## Generating data

check_reliability = function(theta, rtol, atol, max_num_steps) {
  y_hat = solve_sir(theta, rtol, atol, max_num_steps)
  y_hat2 = solve_sir(theta, rtol / 10, atol / 10, max_num_steps)
  max_abs_err = max(abs(y_hat - y_hat2))
  return(max_abs_err)
}

errors = c()
tols = 10^(-c(1:12))
for (tol in tols) {
  errors = c(errors, check_reliability(theta_true, tol, tol, 1e8))
}

qplot(tols, errors, geom = c("point", "line")) +
  scale_x_log10() +
  scale_y_log10()

atol =  1e-6
rtol = 1e-6

dispersion = 5 # noise parameter for negative binomial
mu = solve_sir(theta_true, atol, rtol, 1e8)[, 2]
y = stats::rnbinom(length(ts), mu = mu, size = dispersion)

tibble(t = ts, mu = mu, y = y) %>%
  ggplot() +
  geom_line(aes(t, mu), col = "firebrick") +
  geom_point(aes(t, y)) +
  xlab("Day") +
  ylab("Infected people") +
  ggtitle("Simulated data as points \nUnderlying solution as lines")

rtol_low = 1e-4
atol_low = 1e-3
max_num_steps_low = 100

fit = rstan::sampling(model, list(
  N = length(ts),
  M = M,
  ts = ts,
  y = y,
  initial_conditions = initial_conditions,
  rtol = rtol_low,
  atol = atol_low,
  max_num_steps = max_num_steps_low))

print(get_elapsed_time(fit))
print(fit, pars = c("beta", "gamma", "phi"))


### Tuning the reference method

atol_high = 1e-6
rtol_high = 1e-6
max_num_steps_high = 1e6
draws = rstan::extract(fit, pars = c("beta", "gamma", "phi"))
phi_draws = draws$phi
theta_draws = cbind(draws$beta, draws$gamma)
num_draws = length(phi_draws)

# Compute differences
errors = c()
for (i in 1:num_draws) {
  mae = check_reliability(theta_draws[i,], atol_high, rtol_high,
                          max_num_steps_high)
  errors = c(errors, mae)
}

qplot(errors, geom = "histogram")

### Computing importance weights

log_lh_low = rep(0, num_draws)
log_lh_high = rep(0, num_draws)
for (i in seq_len(num_draws)) {
  y_hat_low = solve_sir(theta_draws[i,], rtol_low, atol_low,
                        max_num_steps_low)
  y_hat_high = solve_sir(theta_draws[i,], rtol_high, atol_high, 
                         max_num_steps_high)
  log_lh_low[i] = sum(dnbinom(y, size = phi_draws[i],
                              mu = y_hat_low[, 2], log = TRUE))
  log_lh_high[i] = sum(dnbinom(y, size = phi_draws[i],
                               mu = y_hat_high[, 2], log = TRUE))
}
log_weights = log_lh_high - log_lh_low

qplot(log_weights, geom = "histogram")

### Computing $\hat{k}$ diagnostic

loo::psis(log_weights)
