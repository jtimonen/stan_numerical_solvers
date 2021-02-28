library(cmdstanr)
library(posterior)
library(loo)
source("functions.R")

# Compile model
model <- cmdstan_model("../stan/lv.stan", include_paths = "../stan")
model_sim <- cmdstan_model("../stan/lv_sim.stan", include_paths = "../stan")

# Load test data
dat <- readRDS(file = "../data/data_lv.rds")

# Fit model
h <- 1.0
tol <- 1e-6
fit <- sample_lv(model, dat, h, tol = tol,
  solver = 0, refresh = 1000, chains = 10, init = 0
)

# Plot solution with posterior mean params
theta_mean <- apply(fit$draws(variables = "theta"), 3, mean)
t_eval <- seq(0.1, 12, by = 0.1)
out <- solve_lv(model_sim, theta_mean, t_eval, h)
plot_lv(dat, out)
diag <- fit$cmdstan_diagnose()

# Comparison
fit_ref <- sample_lv(model, dat, h, tol = tol,
  solver = 1, refresh = 1000, chains = 10, init = 0
)

# Print info
runtime_info(fit)
runtime_info(fit_ref)

# Get log ratios, compute r_eff and pareto_k
#log_ratios <- fit$draws(variables = "log_ratio")
#r_eff <- as.numeric(relative_eff(log_ratios))
#psis <- psis(log_ratios, r_eff = r_eff)

#msg3 <- paste0("Pareto_k = ", round(psis$diagnostics$pareto_k, 3), "\n")
#cat(msg3)

