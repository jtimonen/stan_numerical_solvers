# Create a block of Stan code
stan_block <- function(block_name, block_code) {
  block_name <- trimws(tolower(block_name))
  block_code <- trimws(block_code, whitespace = "[\n]")
  paste0(block_name, " {\n", block_code, "\n}\n")
}

# Create many blocks of Stan code
stan_code <- function(blocks, codes) {
  model_code <- ""
  J <- length(blocks)
  for (j in seq_len(J)) {
    model_code <- paste(model_code, stan_block(blocks[j], codes[j]), sep = "\n")
  }
  return(model_code)
}

# Create a CmdStanModel for prior sampling
create_prior_model <- function(pars, tpars, prior) {
  cat("* Creating CmdStanModel for prior sampling...\n")
  blocks <- c("parameters", "transformed parameters", "model")
  codes <- c(pars, tpars, prior)
  model_code <- stan_code(blocks, codes)
  cat(model_code)
  tmp_file <- cmdstanr::write_stan_file(model_code)
  cmdstanr::cmdstan_model(tmp_file)
}

# Create a CmdStanModel for simulating data given param draws
create_simulator_model <- function(funs, data, tdata, pars, tpars, ode, gq) {
  cat("* Creating CmdStanModel for simulating...\n")
  blocks <- c(
    "functions", "data", "transformed data",
    "parameters", "transformed parameters", "generated quantities"
  )
  gq <- paste(ode, gq, sep="\n")
  codes <- c(funs, data, tdata, pars, tpars, gq)
  model_code <- stan_code(blocks, codes)
  cat(model_code)
  tmp_file <- cmdstanr::write_stan_file(model_code)
  cmdstanr::cmdstan_model(tmp_file)
}

# Create a CmdStanModel for posterior sampling
create_posterior_model <- function(funs, data, tdata, obsdata, 
                                   pars, tpars, prior, ode, lik) {
  cat("* Creating CmdStanModel for posterior sampling...\n")
  blocks <- c(
    "functions", "data", "transformed data",
    "parameters", "transformed parameters", "model"
  )
  post <- paste(prior, ode, lik, sep = "\n")
  data <- paste(data, obsdata, sep="\n")
  codes <- c(funs, data, tdata, pars, tpars, post)
  model_code <- stan_code(blocks, codes)
  cat(model_code)
  tmp_file <- cmdstanr::write_stan_file(model_code)
  cmdstanr::cmdstan_model(tmp_file)
}

# Create all CmdStanModels
create_cmdstan_models <- function(funs, data, tdata, obsdata, pars,
                                  tpars, prior, ode, lik, gq) {
  list(
    prior = create_prior_model(pars, tpars, prior),
    sim = create_simulator_model(funs, data, tdata, pars, tpars, ode, gq),
    post = create_posterior_model(funs, data, tdata, obsdata, 
                                  pars, tpars, prior, ode, lik)
  )
}
