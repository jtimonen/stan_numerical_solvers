require(rstan)

# Create a stan model containing only the functions
two_spaces <- "  " 
f_list <- readLines('functions.stan')
functions <- paste(unlist(f_list), collapse = paste0("\n", two_spaces))
functions <- paste0(two_spaces, functions)
model_code <- paste(c("functions {", functions, "}\n"), collapse = "\n")
model <- rstan::stan_model(model_code = model_code)

# Expose stan functions to R
rstan::expose_stan_functions(model)

# Load solve_tridiag for comparison
source('functions.R')

# Test data
a <- c(1,1,1)
b <- c(2,2,2,2)
d <- runif(4)

# Check if the stan functions work similarly as the R counterpart
x1 <- solve_tridiag(a, b, a, d)
x2 <- stan_solve_tridiag(a, b, a, d)
x3 <- stan_solve_tridiag_sym(a, b, d)
print(x1)
print(x2)
print(x3)

# yay, they do!
