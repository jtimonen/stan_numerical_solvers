require(rstan)

# Compile diffusion.stan
model <- rstan::stan_model('diffusion.stan')

# Expose stan_solve_tridiag to R
rstan::expose_stan_functions(model)

# Load solve_tridiag
source('test_functions.R')

# Test data
a <- c(1,1,1)
b <- c(2,2,2,2)
d <- runif(4)

# Check if they work similarly
x1 <- solve_tridiag(a, b, a, d)
x2 <- stan_solve_tridiag(a, b, a, d)
print(x1)
print(x2)

# yay, they do!
