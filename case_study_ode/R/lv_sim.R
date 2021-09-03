library(cmdstanr)
source("functions.R")

# Compile model
model_sim <- cmdstan_model("../stan/lv_sim.stan", include_paths = "stan")

# Simulate data
N <- 30
D <- 2
theta <- c(1.2, 0.8)
t_data <- seq(12 / N, 12, length.out = N) + 0.3 * rnorm(n = length(N))
if (min(t_data) < 0) {
  t_data <- sort(t_data) - min(t_data) + 0.01
}
h <- 1.0
out <- solve_lv(model_sim, theta, t_data, h)

# Add noise
y_data <- out$y_ref
y_data <- y_data + matrix(0.3 * rnorm(n = D * N), N, D)
dat <- list(t_data = t_data, y_data = y_data, N = N, D = D)

# Plot
plot_lv(dat, out)

# Write to file
saveRDS(file = "../data/data_lv.rds", dat)
