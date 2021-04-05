library(cmdstanr)

# Compile model
model <- cmdstan_model("stan/lv.stan")

# Load test data
dat <- readRDS(file = "test_data.rds")

# Subsample to N data points (out of 3000)
N <- 20
inds <- sort(sample.int(length(dat$t_data), N, replace = FALSE))
t_data <- dat$t_data[inds]
y_data <- dat$y_data[inds, ]

# Fit models
y0 <- c(1, 2)
t0 <- 0
stan_data <- list(
  N = N,
  t_eval = t_data,
  y0 = y0,
  t0 = t0,
  y_data = y_data,
  REL_TOL = 1e-6,
  ABS_TOL = 1e-6,
  MAX_NUM_STEPS = 30
)

fit <- model$sample(data = stan_data, refresh = 1, chains = 4)

# Print results
t <- fit$time()$total
cat("\nN =", N, "\n")
msg <- paste0("Total time: ", round(t, 3), " s\n")
cat(msg)
