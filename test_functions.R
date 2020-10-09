#' Forward Euler method
#' 
#' @param u_init initial state u(0, x)
#' @param dt discretization step in t
#' @param dx discretization step in x
#' @param Nt number of time steps to take
#' @param Kappa diffusion constant
solve_fe <- function(u_init, dt, dx, Nt, Kappa) {
  Nx <- length(u_init)
  U <- matrix(0, Nt + 1, Nx)
  U[1, ] <- u_init
  f <- Kappa * dt / (dx**2)
  for (n in seq_len(Nt)) {
    u_n <- U[n,]
    u_np1 <- u_n
    for (i in 2:(Nx - 1)) {
      u_np1[i] <-  u_n[i] + f * (u_n[i + 1] - 2 * u_n[i] + u_n[i - 1]) 
    }
    U[n + 1,] <- u_np1
    U[n + 1, 0] <- 0
    U[n + 1, Nx] <- 1
  }
  return(U)
}

# Function to plot the solution
plot_u <- function(x, U) {
  lwd <- 2
  plot(x, U[1,], type = 'l', main = 'Solution of u(t,x)',
       ylab = 'u(t,x)', xaxt = "n", lwd = lwd)
  L <- nrow(U)
  for (n in 2:L) {
    lines(x, U[n,], lwd = lwd) 
  }
  #legend(0.8, 0.4, legend = c("t = 0", "t = T"), col = c(col1, col2),
  #       lty = c(1,1), lwd = c(lwd, lwd))
  axis(1, at = c(0, 0.5, 1.0), labels = c("0", "L/2", "L"))
}
