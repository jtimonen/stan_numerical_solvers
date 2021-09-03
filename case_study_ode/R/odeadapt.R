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
  gq <- paste(ode, gq, sep = "\n")
  codes <- c(funs, data, tdata, pars, tpars, gq)
  model_code <- stan_code(blocks, codes)
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
  data <- paste(data, obsdata, sep = "\n")
  codes <- c(funs, data, tdata, pars, tpars, post)
  model_code <- stan_code(blocks, codes)
  tmp_file <- cmdstanr::write_stan_file(model_code)
  cmdstanr::cmdstan_model(tmp_file)
}

# Create all CmdStanModels
create_cmdstan_models <- function(funs, data, tdata, obsdata, pars,
                                  tpars, prior, ode, lik, gq) {
  list(
    prior = create_prior_model(pars, tpars, prior),
    sim = create_simulator_model(funs, data, tdata, pars, tpars, ode, gq),
    post = create_posterior_model(
      funs, data, tdata, obsdata,
      pars, tpars, prior, ode, lik
    )
  )
}

# Function for simulating ODE solutions and data given parameter( draws)s
simulate <- function(model, params, data, solver_args = list()) {
  stopifnot(is(model, "CmdStanModel"))
  stopifnot(is(params, "draws"))
  stopifnot(is(data, "list"))
  stopifnot(is(solver_args, "list"))
  if(is.null(solver_args$RTOL)) solver_args$RTOL <- 1e-6
  if(is.null(solver_args$ATOL)) solver_args$ATOL <- 1e-6
  if(is.null(solver_args$MAX_NUM_STEPS)) solver_args$MAX_NUM_STEPS <- 1e6
  model$generate_quantities(
    data = c(data, solver_args), 
    fitted_params = params
  )
}

# Using simulate with different tolerances
simulate_many <- function(model, params, data, 
                          atols, rtols, MAX_NUM_STEPS=NULL) {
  stopifnot(is(params, "draws"))
  J1 <- length(atols)
  J2 <- length(rtols)
  S <- niterations(params) * nchains(params)
  TIME <- array(0, dim=c(J1, J2))
  XSIM <- array(0, dim=c(J1, J2, S, data$N*data$D))
  for(j1 in 1:J1) {
    for(j2 in 1:J2) {
      solver_args <- list(
        ATOL = atols[j1], 
        RTOL=rtols[j2],
        MAX_NUM_STEPS = MAX_NUM_STEPS
      )
      sim <- simulate(model, params, data, solver_args)
      XSIM[j1, j2,,]  <- merge_chains(sim$draws("x"))[,1,, drop=TRUE]
      TIME[j1, j2] <- sim$time()$total
    }
  }
  return(list(times=TIME,sims=XSIM))
}

# Compute error to most accurate solution
compute_errors <- function(XSIM, fun="max") {
  J1 <- dim(XSIM)[1]
  J2 <- dim(XSIM)[2]
  ERR <- array(0, dim=c(J1,J2))
  I_sim_best <- XSIM[1,1,,]
  for(j1 in 1:J1) {
    for(j2 in 1:J2) {
      I_sim <- XSIM[j1,j2,,]
      abs_errors <- abs(as.vector(I_sim_best)-as.vector(I_sim))
      ERR[j1,j2] <- eval(call(fun, abs_errors))
    }
  }
  return(ERR)
}

# Runtimes plot
plot_sim_times <- function(atols, rtols, TIME) {
  par(mfrow=c(1,2))
  image(log10(atols), log10(rtols), TIME, main = "time (s)")
}


# Errors plot
plot_sim_errors <- function(atols, rtols, ERR) {
  par(mfrow=c(1,2))
  image(log10(atols), log10(rtols), log(ERR), main="log error")
}
