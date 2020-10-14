require(rstan)
require(ggplot2)

model <- rstan::stan_model(file = 'sir.stan')
expose_stan_functions(model)



# Fit model

