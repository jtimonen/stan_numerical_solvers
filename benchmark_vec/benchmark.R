library(cmdstanr)

# Compile models
mod1 <- cmdstan_model("stan/lv_loop.stan")
mod2 <- cmdstan_model("stan/lv_vec.stan")

# Load test data
dat <- readRDS(file = "test_data.rds")

# Subsample to N data points (out of 3000)
N <- 400
inds <- sort(sample.int(length(dat$t_data), N, replace = FALSE))
t_data <- dat$t_data[inds]
y_data <- dat$y_data[inds, ]

# Fit models
y0 <- c(1, 2)
t0 <- 0
stan_data <- list(N = N, t_eval = t_data, y0 = y0, t0 = t0, y_data = y_data)
fit1 <- mod1$sample(data = stan_data, refresh = 0, init = 0)
fit2 <- mod2$sample(data = stan_data, refresh = 0, init = 0)

# Print results
t1 <- fit1$time()$total
t2 <- fit2$time()$total
cat("\nN =", N, "\n")
msg1 <- paste0("Time using looped likelihood: ", round(t1, 3), " s\n")
msg2 <- paste0("Time using vectorized likelihood: ", round(t2, 3), " s\n")
prc <- 100.0 * (t1 - t2) / t1
msg3 <- paste0("Improvement: ", round(prc, 3), "%")
cat(msg1)
cat(msg2)
cat(msg3)
