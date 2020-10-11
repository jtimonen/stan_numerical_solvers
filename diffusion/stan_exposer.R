# Create a stan model containing only the functions
stan_exposer <- function(parent_dir = NULL){
  two_spaces <- "  " 
  STAN_FILES <- c("be_tridiag.stan")
  if (!is.null(parent_dir)) {
    STAN_FILES <- file.path(parent_dir, STAN_FILES)
  }
  sf <- paste(STAN_FILES, collapse = ", ")
  msg <- paste0("Reading '", sf, "'... ")
  cat(msg)
  f_list <- lapply(STAN_FILES, FUN = readLines)
  functions <- paste(unlist(f_list), collapse = paste0("\n", two_spaces))
  functions <- paste0(two_spaces, functions)
  model_code <- paste(c("functions {", functions, "}\n"), collapse = "\n")
  model <- rstan::stan_model(model_code = model_code)
  cat("Done.\n")
  return(model)
}
