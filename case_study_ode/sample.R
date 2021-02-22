library(cmdstanr)
library(posterior)
library(loo)
source("functions.R")

# Compile model
model <- cmdstan_model("lv.stan", include_paths = ".")
model_sim <- cmdstan_model("lv_sim.stan", include_paths = ".")

# Load test data
dat <- readRDS(file = "data_small.rds")

# Fit model
h <- 0.1
fit <- sample_lv(model, dat, h, refresh = 1000, chains = 10,
                 parallel_chains = 4)

# Plot solution with posterior mean params
theta_mean <- apply(fit$draws(variables = "theta"), 3, mean)
t_eval <- seq(0.1, 12, by = 0.1)
out <- solve_lv(model_sim, theta_mean, t_eval, h)
plot_lv(dat, out)
diag <- fit$cmdstan_diagnose()

# Get log ratios, compute r_eff and pareto_k
log_ratios <- fit$draws(variables = "log_ratio")
r_eff <- as.numeric(relative_eff(log_ratios))
psis <- psis(log_ratios, r_eff = r_eff)

# Print info
str1 <- round(fit$time()$total, 3)
str2 <- round(psis$diagnostics$pareto_k, 3)
msg <- paste0("Time = ", str1, " s, pareto_k = ", str2, "\n")
cat(msg)

