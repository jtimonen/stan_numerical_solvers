library(cmdstanr)
source("functions.R")

# Compile model
model <- cmdstan_model("lv.stan", include_paths = ".")
model_sim <- cmdstan_model("lv_sim.stan", include_paths = ".")

# Load test data
dat <- readRDS(file = "data_small.rds")

# Fit model
h <- 1.0
fit <- sample_lv(model, dat, h)

# Plot solution with posterior mean params
theta_mean <- apply(fit$draws(variables = "theta"), 3, mean)
t_eval <- seq(0.1, 12, by = 0.1)
out <- solve_lv(model_sim, theta_mean, t_eval, h)
plot_lv(dat, out)

