#' Solve a tridiagonal linear system Ax = d
#'
#' @param a lower diagonal of A
#' @param b diagonal of matrix A
#' @param c upper diagonal of A
#' @param d right-hand side of the system
solve_tridiag <- function(a, b, c, d) {
  n <- length(b)
  x <- rep(0, n)
  
  # Forward sweep
  for (i in 2:n) {
    w <- a[i - 1] / b[i - 1]
    b[i] <- b[i] - w*c[i - 1]
    d[i] <- d[i] - w*d[i - 1]
  }
  
  # Back substitution
  x[n] <- d[n]/b[n]
  for (i in (n - 1):1) {
    x[i] <- (d[i] - c[i]*x[i + 1]) / b[i]
  }
  return(x)
}

#' Solve a symmetric tridiagonal linear system Ax = d
#'
#' @param a upper/lower diagonal of A
#' @param b diagonal of matrix A
#' @param d right-hand side of the system
solve_tridiag_sym <- function(a, b, d) {
  n <- length(b)
  x <- rep(0, n)
  
  # Forward sweep
  for (i in 2:n) {
    w <- a[i - 1] / b[i - 1]
    b[i] <- b[i] - w*a[i - 1]
    d[i] <- d[i] - w*d[i - 1]
  }
  
  # Back substitution
  x[n] <- d[n]/b[n]
  for (i in (n - 1):1) {
    x[i] <- (d[i] - a[i]*x[i + 1]) / b[i]
  }
  return(x)
}

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
  
  # Create the diagonal, and upper and lower secondary diagonals of A
  A_diag <- rep(1.0 + 2.0*K, Nx)
  A_diag[1] <- 1.0 + K
  A_diag[Nx] <- 1.0 + K
  A_band <- rep(-K, Nx - 1)
  
  for (n in seq_len(Nt)) {
    u_n <- U[n,]
    #u_np1 <- solve_tridiag(A_band, A_diag, A_band, u_n)
    u_np1 <- solve_tridiag_sym(A_band, A_diag, u_n)
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
#' @param main main title
plot_u <- function(x, U, T_max, main) {
  lwd <- 2
  col1 <- "gray30"
  col2 <- "firebrick3"
  col3 <- "orange"
  leg <- c("t = 0", paste("t =", T_max/4), paste("t =", T_max))
  plot(x, U[1,], col = col1, type = 'l', main = main,
       ylab = 'u(t,x)', xaxt = "n", lwd = lwd)
  L <- nrow(U)
  idx <- round(L/4)
  lines(x, U[idx,], lwd = lwd, col = col2)
  lines(x, U[L,], lwd = lwd, col = col3) 
  legend(0.7, 0.4, legend = leg, 
         col = c(col1, col2, col3), lty = c(1,1,1), lwd = c(lwd, lwd, lwd))
  axis(1, at = c(0, 0.5, 1.0), labels = c("0", "L/2", "L"))
}


#' Function to plot the development of the solution
#' 
#' @param U matrix of solutions
plot_U <- function(U, T_max) {
  image(U, main = 'Solution of u(t,x)', xlab = 't', ylab = 'x', xaxt = "n")
  axis(1, at = c(0, 1.0), labels = c("0", T_max))
}
