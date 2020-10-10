require(rstan)

source('test_functions.R')

# Setup
L <- 1
T_max <- 0.1  # time interval length
Kappa <- 0.2  # true value of kappa
dx <- 1e-2    # discretization step in x
dt <- 1e-4    # discretization step in t

# Define initial heat distribution u(t=0, x)
x_grid <- seq(0, 1, by = dx)
u_init <- as.numeric(x_grid > (0.5 * L))

# Solve u(t=T, x) and plot
Nt <- T_max/dt
U <- solve_be(u_init, dt, dx, Nt, Kappa)
main <- paste0('Solution with Kappa=', Kappa, ', dx=', dx, ", dt=", dt)
plot_u(x_grid, U, T_max, main)
