
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
solve_lv <- function(model, theta, t_eval, h) {
  
  # Create input
  y0 <- c(1, 2)
  t0 <- 0
  N <- length(t_eval)
  num_steps <- ceiling(max(t_eval) / h)
  stan_data <- list(
    N = N, t_eval = t_eval, y0 = y0, t0 = t0,
    h = h, num_steps = num_steps
  )
  
  # Run Stan
  fit <- model$sample(
    data = stan_data,
    fixed_param = TRUE,
    iter_sampling = 1,
    chains = 1,
    init = list(list(theta = theta)),
    refresh = 0
  )
  
  # Get output
  list(
    y_grid_rk4 = get_output(fit, "y_grid_rk4", 1, 1),
    y_rk4 = get_output(fit, "y_rk4", 1, 1),
    y_ref = get_output(fit, "y_ref", 1, 1),
    theta = as.vector(fit$draws(variables = "theta")[1,1,]),
    t_grid = seq(0, num_steps * h, by = h)
  )
}

# Function
sample_lv <- function(model, data, h, ...){
  
  # Create input
  y0 <- c(1, 2)
  t0 <- 0
  t_data <- data$t_data
  y_data <- data$y_data
  N <- length(t_data)
  num_steps <- ceiling(max(t_data) / h)
  stan_data <- list(
    N = N, t_eval = t_data, y0 = y0, t0 = t0,
    h = h, num_steps = num_steps, y_data = y_data
  )
  
  # Run Stan
  fit <- model$sample(data = stan_data, ...)
  return(fit)
}

# Function
plot_lv <- function(dat, out) {
  
  par(mfrow = c(2, 1))
  par(mar = c(2.25, 4, 1, 1))
  cols <- c("steelblue", "firebrick3")
  T_max <- max(dat$t_data)
  for (j in 1:2) {
    plot(dat$t_data, dat$y_data[, j],
         col = "black", pch = 16,
         xlim = c(0, T_max), ylim = c(0, 3.6),
         xlab = "t", ylab = paste0("y", j)
    )
    lines(t_eval, out$y_ref[, j], col = cols[1])
    lines(t_eval, out$y_rk4[, j], col = cols[2], lty = 2)
    points(out$t_grid, out$y_grid_rk4[, j], col = cols[2], pch = 4)
    legend(
      x = 9.5, y = 3.5, lty = c(1, 2), lwd = 2, col = cols,
      legend = c("rk45", "rk4")
    )
  }
}
