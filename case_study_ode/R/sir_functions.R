# Helper function
solve_sir <- function(model, t_eval, rtol, atol, max_num_steps, theta){
  stan_data <- list(N = length(t_eval), 
                    t_eval = t_eval, 
                    RTOL = rtol, 
                    ATOL = atol, 
                    MAX_NUM_STEPS = max_num_steps, 
                    S = nrow(theta), 
                    THETA = theta,
                    pop_size = 1000,
                    I0 = 10
  )
  print(stan_data)
  out <- model$sample(data = stan_data, fixed_param = TRUE, iter_sampling = 1)
  return(out)
}
