library(cmdstanr)
source("sir_functions.R")

# Compile model
fp <-  "../stan/sir_fixedparam.stan"
model <- cmdstan_model(fp, include_paths = "../stan")

# Simulate data
t_data <- seq(1, 40, by = 1)
rtol <- 1e-5
atol <- 1e-5
max_num_steps <- 1e6
THETA <- matrix(c(1,1,0.5,1,0.7,1), 3, 2)

# ASd
fit <- solve_sir(model, t_data, rtol, atol, max_num_steps, THETA)
print(dim(fit$draws(variables = "y_hat")))
