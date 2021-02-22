library(cmdstanr)
source("functions.R")

# Compile model
model_sim <- cmdstan_model("lv_sim.stan",
  include_paths = ".",
  force_recompile = F
)

# Simulate data
N <- 10
D <- 2
theta <- c(1.2, 0.8)
t_data <- seq(1, 12, length.out = N) + 0.3*rnorm(n = length(N))
h <- 0.1
out <- solve_lv(model_sim, theta, t_data, h)

# Add noise
y_data <- out$y_ref
y_data <- y_data + matrix(0.3*rnorm(n = D*N), N, D)
dat <- list(t_data = t_data, y_data = y_data, N = N, D = D)

# Plot
plot_lv(dat, out)

# Write to file
saveRDS(file = "data_small.rds", dat)
