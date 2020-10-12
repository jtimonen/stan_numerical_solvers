require(rstan)
require(ggplot2)
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html


model <- rstan::stan_model(file = 'sir.stan')
expose_stan_functions(model)


# Setup
N <- 1000
I0 <- 10
y0 <- c(N - I0, I0, 0)
t0 <- 0
ts <- seq(1, 30, by = 1)
theta <- c(1,0.2)
rtol <- 1e-10
atol <- 1e-10
max_steps <- 1e6

# Solve ODE
odesolve <- function(y0, t0, ts, theta, N, rtol, atol, max_steps){
  x_r <- numeric(0)
  y_hat <- stan_odesolve(y0, t0, ts, theta, x_r, N, rtol, atol, max_steps)
  L <- length(y_hat)
  y_hat <- unlist(y_hat)
  n <- length(y_hat)
  y_hat <- matrix(y_hat, L, n/L, byrow = TRUE)
  return(y_hat)
}

y_hat <- odesolve(y0, t0, ts, theta, N, rtol, atol, max_steps)
yyy <- rbind(y0, y_hat)
ttt <- c(t0, ts)
# Plot
y <- as.vector(yyy)
day <- rep(ttt, 3)
comp <- rep( c("S", "I", "R"), each=length(ttt))
df <- data.frame(day, y, as.factor(comp))
colnames(df) <- c("Day", "y", "Compartment")

aes <- aes_string(x = "Day", y = "y", group = "Compartment", color = "Compartment")
plt <- ggplot(df, aes) + geom_line(lwd = 1) + geom_point()

