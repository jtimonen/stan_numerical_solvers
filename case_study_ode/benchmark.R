library(cmdstanr)
library(posterior)
library(loo)
source("functions.R")

# Compile models
mod1 <- cmdstan_model("stan/lv_doc.stan", include_paths = "stan")
mod2 <- cmdstan_model("stan/lv_vec.stan", include_paths = "stan")

# Load test data (created using sim.py)
dat <- readRDS(file = "data_small.rds")

# Fit model
y0 <- c(1, 2)
t0 <- 0
t_data <- dat$t_data
y_data <- dat$y_data
N <- length(t_data)

# Fit first model
stan_data <- list(N = N, t_eval = t_data, y0 = y0, t0 = t0, y_data = y_data)
fit1 <- mod1$sample(data = stan_data, refresh = 1000, init = 0)
fit2 <- mod2$sample(data = stan_data, refresh = 1000, init = 0)

t1 <- fit1$time()$total
t2 <- fit2$time()$total

cat(paste0(t1, "\n", t2, "\n"))
