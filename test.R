require(rstan)

source('test_functions.R')

# Setup
L <- 1        # rod length
T_end <- 1.0  # time interval length
Kappa <- 0.25 # true value of kappa
dx <- 1e-3    # discretization step in x
dt <- 1e-4    # discretization step in t

# Define initial heat distribution u(t=0, x)
x_grid <- seq(0, 1, by = dx)
u_init <- x_grid > (0.5 * L)

# Solve u(t=T, x) and plot
u_end <- solve_u(u_init, T_end, dt, dx)
plot_u(x_grid, u_init, u_end)
