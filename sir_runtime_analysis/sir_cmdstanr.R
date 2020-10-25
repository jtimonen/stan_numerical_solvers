# Requirements
require(cmdstanr)
library(tidyverse)
require(bayesplot)
require(ggplot2)
require(loo)
require(stats)
require(posterior)

# Build model and expose functions
fp <- file.path("..", "stan", "sir.stan")
mod <- cmdstan_model(stan_file = fp)

# Define data
N <- 16                                # number of days measured
M <- 1000                              # population size
I0 <- 5                                # number of infected on day 0
y0 <- c(M - I0, I0, 0)                 # S, I, R on day 0
ts <- seq_len(N)                       # measurement times
y <- stats::dlnorm(ts, meanlog = 2, sdlog = 0.5) * 4000
y <- round(y)

atol <- 1e-5
rtol <- 1e-5
max_num_steps <- 1e8

dat = list(
  N = length(ts),
  M = M,
  ts = ts,
  y = y,
  initial_conditions = y0,
  rtol = atol,
  atol = atol,
  max_num_steps = max_num_steps
)

# Fit model
a <- mod$sample(data = dat)
