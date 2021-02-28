library(cmdstanr)
library(ggplot2)
source("sir_functions.R")

# Compile model
fp <-  "../stan/sir_fixedparam.stan"
model_fp <- cmdstan_model(fp, include_paths = "../stan")

# Setup
t_data <- seq(1, 40, by = 1)
rtol <- 1e-8
atol <- 1e-8
max_num_steps <- 1e6
THETA <- matrix(c(0.5, 0.2), 1, 2, byrow = TRUE)

# Solve ODE
df <- solve_sir(model_fp, t_data, rtol, atol, max_num_steps, THETA)

# Plot
plt <- ggplot(df, aes(x = t, y = y_hat, group = draw, color = draw)) +
  geom_line() + facet_grid(.~ var)

# Get only infected
dat <- df[df$var == "I",]
t_data <- dat$t
x_data <- dat$y_hat

# Add NB noise
phi <- 5
y_data <- rnbinom(n = length(t_data), size = phi, mu = x_data)

# Save
dat <- list(t_data = t_data, y_data = y_data)
#saveRDS(file = "../data/data_sir.rds", dat)
