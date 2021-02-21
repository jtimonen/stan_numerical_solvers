library(cmdstanr)
source("functions.R")

# Compile model
model_sim <- cmdstan_model("lv_sim.stan",
  include_paths = ".",
  force_recompile = F
)

# Load test data
dat <- readRDS(file = "data_small.rds")

# Plot one draw
theta <- c(1.0, 0.8)
t_eval <- seq(0.1, 12, by = 0.1)
out <- solve_lv(model_sim, theta, t_eval, h)
plot_lv(dat, out)
