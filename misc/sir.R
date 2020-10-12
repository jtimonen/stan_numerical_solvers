
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

require(rstan)

model <- rstan::stan_model(file = 'sir.stan')
expose_stan_functions(model)

N <- 100
f <- stan_sir(0.0, c(1,2,3), c(1,1), 0.0, 10.0)
