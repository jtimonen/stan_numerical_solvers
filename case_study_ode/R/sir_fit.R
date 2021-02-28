library(cmdstanr)
library(posterior)
library(loo)
source("sir_functions.R")

# Compile models
f1 <- "../stan/sir_fixedparam.stan"
f2 <- "../stan/sir.stan"
model_fp <- cmdstan_model(f1, include_paths = "../stan")
model <- cmdstan_model(f2, include_paths = "../stan")

# Load test data
dat <- readRDS(file = "../data/data_sir.rds")
t_data <- dat$t_data
y_data <- dat$y_data

# Setup
rtol <- 1e-6 # 1e-3  much slower
atol <- 1e-6
max_num_steps <- 1e6

# Fit model
N <- length(t_data)
stan_data <- list(N = N, 
                  t_data = t_data, 
                  y_data = y_data,
                  RTOL = rtol,
                  ATOL = atol, 
                  MAX_NUM_STEPS = max_num_steps, 
                  pop_size = 1000,
                  I0 = 15
)
fit <- model$sample(data = stan_data)
t1 <- fit$time()$total
print(t1)

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

