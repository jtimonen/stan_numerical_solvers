library(cmdstanr)

# Compile model
model <- cmdstan_model("lv_simulate.stan",
  include_paths = ".",
  force_recompile = T
)

# Create test input
t_data <- seq(0.1, 16, by = 0.1)
N <- length(t_data)

# Set initial state
y0 <- c(1, 2)
t0 <- 0
h <- 1.0
num_steps <- ceiling(max(t_data) / h)
stan_data <- list(
  N = N, t_data = t_data, y0 = y0, t0 = t0,
  h = h, num_steps = num_steps
)

# Function
get_output <- function(fit, name, chain_idx, draw_idx) {
  y_out <- fit$draws(variables = name)
  yj <- as.matrix(y_out[draw_idx, chain_idx, ])
  N_grid <- length(yj) / 2
  y1 <- yj[1:N_grid]
  y2 <- yj[(N_grid + 1):(2 * N_grid)]
  return(cbind(y1, y2))
}

# Function
fit_model <- function(model, theta) {
  fit <- model$sample(
    data = stan_data,
    fixed_param = TRUE,
    iter_sampling = 1,
    chains = 1,
    init = list(list(theta = theta))
  )
  list(
    y_grid_rk4 = get_output(fit, "y_grid_rk4", 1, 1),
    y_rk4 = get_output(fit, "y_rk4", 1, 1),
    y_ref = get_output(fit, "y_ref", 1, 1),
    theta = as.vector(fit$draws(variables = "theta")[1,1,])
  )
}

# Plot one draw
t_grid <- seq(0, num_steps * h, by = h)
j <- 1

out <- fit_model(model, c(1, 1))

par(mfrow = c(2, 1))
par(mar = c(2.25, 4, 1, 1))
cols <- c("steelblue", "firebrick3")
for (j in 1:2) {
  plot(t_data, out$y_ref[, j], "n",
    col = "gray", pch = 20,
    xlim = c(0, max(t_data)), ylim = c(0, 3.6),
    xlab = "t", ylab = paste0("y", j)
  )
  lines(t_data, out$y_ref[, j], col = cols[1])
  lines(t_data, out$y_rk4[, j], col = cols[2], lty = 2)
  points(t_grid, out$y_grid_rk4[, j], col = cols[2], pch = 20)
  if (j == 1) {
    legend(
      x = 0.5, y = 3.5, lty = c(1, 2), lwd = 2, col = cols,
      legend = c("rk45", "rk4")
    )
  }
}
