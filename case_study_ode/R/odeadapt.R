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

# Create full Stan code for prior sampling model
prior_model_code <- function(pars, tpars, prior) {
  cat("* Creating CmdStanModel for prior sampling...\n")
  blocks <- c("parameters", "transformed parameters", "model")
  codes <- c(pars, tpars, prior)
  stan_code(blocks, codes)
}

# Create full Stan code for simulator model
simulator_model_code <- function(funs, data, tdata, obsdata,
                                 pars, tpars, ode, gq) {
  cat("* Creating CmdStanModel for simulating...\n")
  blocks <- c(
    "functions", "data", "transformed data",
    "parameters", "transformed parameters", "generated quantities"
  )
  data <- paste(data, obsdata, sep = "\n")
  gq <- paste(ode, gq, sep = "\n")
  codes <- c(funs, data, tdata, pars, tpars, gq)
  stan_code(blocks, codes)
}

# Create full Stan code for posterior sampling model
posterior_model_code <- function(funs, data, tdata, obsdata,
                                 pars, tpars, prior, ode, lik) {
  cat("* Creating CmdStanModel for posterior sampling...\n")
  blocks <- c(
    "functions", "data", "transformed data",
    "parameters", "transformed parameters", "model"
  )
  post <- paste(prior, ode, lik, sep = "\n")
  data <- paste(data, obsdata, sep = "\n")
  codes <- c(funs, data, tdata, pars, tpars, post)
  stan_code(blocks, codes)
}

# Create all CmdStanModels
create_cmdstan_models <- function(funs, data, tdata, obsdata, pars,
                                  tpars, prior, ode, lik, gq) {
  codes <- list(
    prior = prior_model_code(pars, tpars, prior),
    simulator = simulator_model_code(
      funs, data, tdata, obsdata,
      pars, tpars, ode, gq
    ),
    posterior = posterior_model_code(
      funs, data, tdata, obsdata,
      pars, tpars, prior, ode, lik
    )
  )
  j <- 0
  models <- list()
  for (code in codes) {
    j <- j + 1
    models[[j]] <- cmdstan_model(write_stan_file(code))
  }
  names(models) <- names(codes)
  return(models)
}

# Function for simulating ODE solutions and data given parameter( draws)s
simulate <- function(model, params, data, solver_args = list()) {
  stopifnot(is(model, "CmdStanModel"))
  stopifnot(is(params, "draws"))
  stopifnot(is(data, "list"))
  stopifnot(is(solver_args, "list"))
  if (is.null(solver_args$RTOL)) solver_args$RTOL <- 1e-6
  if (is.null(solver_args$ATOL)) solver_args$ATOL <- 1e-6
  if (is.null(solver_args$MAX_NUM_STEPS)) solver_args$MAX_NUM_STEPS <- 1e6
  model$generate_quantities(
    data = c(data, solver_args),
    fitted_params = params
  )
}

# Using simulate with different tolerances
simulate_many <- function(model, params, data,
                          atol, rtol, MAX_NUM_STEPS = NULL) {
  stopifnot(is(params, "draws"))
  J1 <- length(atol)
  J2 <- length(rtol)
  S <- niterations(params) * nchains(params)
  TIME <- array(0, dim = c(J1, J2))
  XSIM <- array(0, dim = c(J1, J2, S, data$N * data$D))
  LL <- array(0, dim = c(J1, J2, S))
  for (j1 in 1:J1) {
    for (j2 in 1:J2) {
      solver_args <- list(
        ATOL = atol[j1],
        RTOL = rtol[j2],
        MAX_NUM_STEPS = MAX_NUM_STEPS
      )
      sim <- simulate(model, params, data, solver_args)
      XSIM[j1, j2, , ] <- merge_chains(sim$draws("x"))[, 1, , drop = TRUE]
      TIME[j1, j2] <- sim$time()$total
      LL[j1, j2, ] <- merge_chains(sim$draws("log_lik"))[, 1, 1, drop = TRUE]
    }
  }
  return(list(times = TIME, sims = XSIM, log_liks = LL))
}

# Compute error to most accurate solution
compute_sol_errors <- function(XSIM, fun = "max") {
  J1 <- dim(XSIM)[1]
  J2 <- dim(XSIM)[2]
  ERR <- array(0, dim = c(J1, J2))
  I_sim_best <- XSIM[1, 1, , ]
  for (j1 in 1:J1) {
    for (j2 in 1:J2) {
      I_sim <- XSIM[j1, j2, , ]
      abs_errors <- abs(as.vector(I_sim_best) - as.vector(I_sim))
      ERR[j1, j2] <- eval(call(fun, abs_errors))
    }
  }
  return(ERR)
}

# Compute error in likelihood, compared to most accurate solution
compute_loglik_errors <- function(LL, fun = "max") {
  J1 <- dim(LL)[1]
  J2 <- dim(LL)[2]
  ERR <- array(0, dim = c(J1, J2))
  loglik_best <- LL[1, 1, ]
  for (j1 in 1:J1) {
    for (j2 in 1:J2) {
      loglik <- LL[j1, j2, ]
      abs_errors <- abs(as.vector(loglik_best) - as.vector(loglik))
      ERR[j1, j2] <- eval(call(fun, abs_errors))
    }
  }
  return(ERR)
}

# Runtimes plot
plot_sim_times <- function(atol, rtol, TIME) {
  par(mfrow = c(1, 2))
  plot(log10(rtol), diag(TIME),
    xlab = "log10(tol)", ylab = "time (s)",
    type = "o", pch = 16
  )
  grid()
  image(log10(atol), log10(rtol), TIME, main = "time (s)")
}


# Errors plot
plot_sim_errors <- function(atol, rtol, ERR, log = TRUE) {
  main <- deparse(substitute(ERR))
  if (log) {
    ERR <- log(ERR)
    main <- paste0("log(", main, ")")
  }
  par(mfrow = c(1, 2))
  plot(log10(rtol), diag(ERR),
    xlab = "log10(tol)", ylab = main,
    type = "o", pch = 16
  )
  grid()
  image(log10(atol), log10(rtol), ERR, main = main)
}
