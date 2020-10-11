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


#' Backward Euler method
#' 
#' @param u_init initial state u(0, x)
#' @param dt discretization step in t
#' @param dx discretization step in x
#' @param Nt number of time steps to take
#' @param K diffusion constant
#' @param ul value at left boundary
#' @param ur value at right boundary
be <- function(u_init, dt, dx, Nt, K, ul, ur) {
  Nx <- length(u_init)
  U <- matrix(0, Nt + 1, Nx)
  U[1, ] <- u_init
  K_star <- K * dt / (dx**2)
  msg <- paste0("K_star = ", K_star, "\n")
  cat(msg)
  
  # Create the diagonal, and upper and lower secondary diagonals of A
  A_diag <- rep(1.0 + 2.0*K_star, Nx)
  A_band <- rep(-K_star, Nx - 1)
  
  for (n in seq_len(Nt)) {
    b <- U[n,]
    #u_np1 <- solve_tridiag(A_band, A_diag, A_band, u_n)
    b[1] = b[1] + ul * K_star
    b[Nx] = b[Nx] + ur * K_star
    U[n + 1,] <- solve_tridiag_sym(A_band, A_diag, b)
  }
  return(U)
}

#' Function to plot the solution(s)
#' 
#' @param x grid in x
#' @param U matrix where each row is one solution (R rows)
#' @param t vector where each value is the correspondinng time (length R)
#' @param cols vector of colors (length R)
#' @param main main title
plot_u <- function(x, U, t, main) {
  lwd <- 2
  leg <- paste0("t = ", t)
  plot(x, U[1,], type = 'l', main = main,
       ylab = 'u(t,x)', xaxt = "n", lwd = lwd)
  R <- nrow(U)
  for (r in 2:R) {
    lines(x, U[r,], lwd = lwd)
  }
  legend(0.7, 0.4, legend = leg, lty = rep(1, R), lwd = rep(lwd, R))
  axis(1, at = c(0, 0.5, 1.0), labels = c("0", "L/2", "L"))
}

#' Function to plot the development of the solution
#' 
#' @param U matrix of solutions
plot_U <- function(U, T_max) {
  image(U, main = 'Solution of u(t,x)', xlab = 't', ylab = 'x', xaxt = "n")
  axis(1, at = c(0, 1.0), labels = c("0", T_max))
}
