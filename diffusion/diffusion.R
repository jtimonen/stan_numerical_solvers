source('diffusion/functions.R')

# Setup
L <- 1
T_max <- 1e-1  # time interval length
Kappa <- 0.2  # true value of kappa
dx <- 1e-2    # discretization step in x
dt <- 1e-1    # discretization step in t

# Define initial heat distribution u(t=0, x)
x_grid <- seq(dx, 1 - dx, by = dx)
u_init <- as.numeric(x_grid > (0.5 * L))

# Solve u(t=T, x) and plot
Nt <- T_max/dt
U <- be(u_init, dt, dx, Nt, Kappa, 0.0, 1.0)
main <- paste0('Solution with Kappa=', Kappa, ', dx=', dx, ", dt=", dt)
plot_u(x_grid, U, T_max, main)
