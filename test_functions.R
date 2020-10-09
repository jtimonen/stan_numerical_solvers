# Function to solve u(t,x)
solve_u <- function(u_init, T_end, dt, dx) {
  u <- u_init
  for (t in seq(0, T_end, by = dt)) {
    u_new <- u
    for (i in 2:(length(u) - 1)) {
      f <- (u[i + 1] - 2 * u[i] + u[i - 1]) / (2 * dx)
      u_new[i] <-  u[i] + dt * f
    }
    u <- u_new
  }
  return(u)
}

# Function to plot the solution
plot_u <- function(x, u_init, u_end, col1 = "gray30", col2 = "#198bff") {
  lwd <- 2
  plot(x, u_init, type = 'l', col = col1, main = 'Heat distribution on [0,L]',
       ylab = 'u(t,x)', xaxt = "n", lwd = lwd)
  lines(x, u_end, col = col2, lwd = lwd)
  legend(0.8, 0.4, legend = c("t = 0", "t = T"), col = c(col1, col2),
         lty = c(1,1), lwd = c(lwd, lwd))
  axis(1, at = c(0, 0.5, 1.0), labels = c("0", "L/2", "L"))
}
