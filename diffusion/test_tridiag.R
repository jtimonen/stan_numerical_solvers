require(rstan)
source('functions.R')
source('stan_exposer.R')

# Expose stan functions to R
m <- stan_exposer()
rstan::expose_stan_functions(m)

# Test data
a <- c(1,1,1)
b <- c(2,3,3,2)
d <- runif(4)

# Check if the stan functions work similarly as the R counterpart
x1 <- solve_tridiag(a, b, a, d)
x2 <- stan_solve_tridiag(a, b, a, d)
x3 <- stan_solve_tridiag_sym(a, b, d)
x4 <- stan_solve_tridiag_be(a[1], b, d)
print(x1)
print(x2)
print(x3)
print(x4)

# yay, they do!
