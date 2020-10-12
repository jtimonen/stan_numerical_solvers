require(rstan)

model <- rstan::stan_model(file = 'sir.stan')
expose_stan_functions(model)