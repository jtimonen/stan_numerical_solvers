# Create a stan model containing only the functions
stan_exposer <- function(){
  two_spaces <- "  " 
  STAN_FILES <- c("tridiag.stan", "be.stan")
  f_list <- lapply(STAN_FILES, FUN = readLines)
  functions <- paste(unlist(f_list), collapse = paste0("\n", two_spaces))
  functions <- paste0(two_spaces, functions)
  model_code <- paste(c("functions {", functions, "}\n"), collapse = "\n")
  model <- rstan::stan_model(model_code = model_code)
  return(model)
}
