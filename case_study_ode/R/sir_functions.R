# Helper function
solve_sir <- function(model, t_eval, rtol, atol, max_num_steps, theta) {
  N <- length(t_eval)
  S <- nrow(theta)
  stan_data <- list(
    N = N,
    t_eval = t_eval,
    RTOL = rtol,
    ATOL = atol,
    MAX_NUM_STEPS = max_num_steps,
    S = S,
    THETA = theta,
    pop_size = 1000,
    I0 = 15
  )
  out <- model$sample(data = stan_data, fixed_param = TRUE, iter_sampling = 1)
  draws <- out$draws(variables = "y_hat")
  vars <- c("S", "I", "R")
  D <- length(vars)
  var <- rep(vars, each = N * S)
  time <- rep(rep(t_eval, each = S), D)
  draw <- rep(1:S, N * D)
  df <- data.frame(as.vector(draws), time, as.factor(draw), as.factor(var))
  colnames(df) <- c("y_hat", "t", "draw", "var")
  return(df)
}
