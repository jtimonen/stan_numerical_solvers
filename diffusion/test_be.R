require(rstan)
source('functions.R')
source('stan_exposer.R')

# Expose stan functions to R
m <- stan_exposer()
rstan::expose_stan_functions(m)

# Setup a test case
L <- 1
T_max <- 0.2     # time interval length
K     <- 0.1     # true value of K
dx    <- 1e-2    # discretization step in x
dt    <- 1e-4    # discretization step in t

# Define initial heat distribution u(t=0, x)
x_grid <- seq(0, 1, by = dx)
u_init <- as.numeric(x_grid > (0.5 * L))

# Solve u using the stan and R functions
Nt <- T_max/dt
U <- solve_be(u_init, dt, dx, Nt, K)

u1 <- U[Nt,]
u2 <- stan_be(u_init, dt, dx, T_max, K)

# See if they give same solution
plot(x_grid, u1, ylab = 'u')
lines(x_grid, u2, col = "red")
