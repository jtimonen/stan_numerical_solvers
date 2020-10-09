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
  K <- Kappa * dt / (dx**2)
  msg <- paste0("K = ", K, "\n")
  cat(msg)
  for (n in seq_len(Nt)) {
    u_n <- U[n,]
    u_np1 <- u_n
    for (i in 2:(Nx - 1)) {
      u_np1[i] <-  u_n[i] + K * (u_n[i + 1] - 2 * u_n[i] + u_n[i - 1]) 
    }
    U[n + 1,] <- u_np1
    U[n + 1, 0] <- 0
    U[n + 1, Nx] <- 1
  }
  return(U)
}

#' Create the tridiagonal matrix A
#' 
#' @param Nx number of discretization points in x
#' @param K dimensionless diffusion constant
#' @examples
create_A <- function(Nx, K) {
  A <- matrix(0, Nx, Nx)
  for (i in 2:(Nx - 1)) {
    A[i, i + 1] <- -K
    A[i, i - 1] <- -K
    A[i, i]     <- 1 + 2*K
  }
  A[1, 1]       <- 1.0 + K
  A[Nx, Nx]     <- 1.0 + K
  A[1, 2]       <- -K
  A[Nx, Nx - 1] <- -K
  return(A)
}

#' Backward Euler method
#' 
#' @param u_init initial state u(0, x)
#' @param dt discretization step in t
#' @param dx discretization step in x
#' @param Nt number of time steps to take
#' @param Kappa diffusion constant
solve_be <- function(u_init, dt, dx, Nt, Kappa) {
  Nx <- length(u_init)
  U <- matrix(0, Nt + 1, Nx)
  U[1, ] <- u_init
  K <- Kappa * dt / (dx**2)
  msg <- paste0("K = ", K, "\n")
  cat(msg)
  A <- create_A(Nx, K)
  for (n in seq_len(Nt)) {
    u_n <- U[n,]
    u_np1 <- solve(A, u_n)
    U[n + 1,] <- u_np1
    U[n + 1, 1] <- 0
    U[n + 1, Nx] <- 1
  }
  return(U)
}

#' Function to plot the solution
#' 
#' @param x grid in x
#' @param U matrix of solutions
#' @param T_max max time
plot_u <- function(x, U, T_max) {
  lwd <- 2
  col1 <- "gray30"
  col2 <- "firebrick3"
  col3 <- "orange"
  leg <- c("t = 0", paste("t =", T_max/4), paste("t =", T_max))
  plot(x, U[1,], col = col1, type = 'l', main = 'Solution of u(t,x)',
       ylab = 'u(t,x)', xaxt = "n", lwd = lwd)
  L <- nrow(U)
  idx <- round(L/4)
  lines(x, U[idx,], lwd = lwd, col = col2)
  lines(x, U[L,], lwd = lwd, col = col3) 
  legend(0.7, 0.4, legend = leg, 
         col = c(col1, col2, col3), lty = c(1,1,1), lwd = c(lwd, lwd, lwd))
  #axis(1, at = c(0, 0.5, 1.0), labels = c("0", "L/2", "L"))
}


#' Function to plot the development of the solution
#' 
#' @param U matrix of solutions
plot_U <- function(U, T_max) {
  image(U, main = 'Solution of u(t,x)', xlab = 't', ylab = 'x', xaxt = "n")
  axis(1, at = c(0, 1.0), labels = c("0", T_max))
}
