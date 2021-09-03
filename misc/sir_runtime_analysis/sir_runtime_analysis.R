# Requirements
require(rstan)
library(tidyverse)
require(bayesplot)
require(ggplot2)
require(loo)
require(stats)
require(posterior)

# -------------------------------------------------------------------------

# Build model and expose functions
model <- stan_model("../stan/sir.stan")
expose_stan_functions(model)

N <- 16 # number of days measured
M <- 1000 # population size
I0 <- 20 # number of infected on day 0
initial_conditions <- c(M - I0, I0, 0) # S, I, R on day 0
ts <- seq_len(N) # measurement times
theta_true <- c(1, 0.2) # true parameter values

# Solve the SIR system and format result array
# - theta = [beta, gamma]
solve_sir <- function(theta, rtol, atol, max_num_steps) {
  stan_solve_sir(
    initial_conditions, ts, theta,
    c(0.0), M, rtol, atol, max_num_steps
  ) %>%
    unlist() %>%
    matrix(ncol = 3, byrow = TRUE)
}

# Generate data
atol <- 1e-6
rtol <- 1e-6

dispersion <- 5 # noise parameter for negative binomial
mu <- solve_sir(theta_true, atol, rtol, 1e8)[, 2]
y <- stats::rnbinom(length(ts), mu = mu, size = dispersion)

tibble(t = ts, mu = mu, y = y) %>%
  ggplot() +
  geom_line(aes(t, mu), col = "firebrick") +
  geom_point(aes(t, y)) +
  xlab("Day") +
  ylab("Infected people") +
  ggtitle("Simulated data as points \nUnderlying solution as lines")

# -------------------------------------------------------------------------


TOL <- 10^(-3:-5) # different values to use for atol and rtol
NTOL <- length(TOL)
MNS <- 10^(2:6) # different values to use for max_num_steps
L <- length(MNS)
C <- 10 # number of chains
runtimes <- array(0, c(NTOL, L, C)) # store runtimes here

for (i in seq_len(NTOL)) {
  for (j in seq_len(L)) {
    msg <- paste0("atol=rtol=", TOL[i], ", max_num_steps=", MNS[j], "\n")
    cat(msg)

    dat <- list(
      N = length(ts),
      M = M,
      ts = ts,
      y = y,
      initial_conditions = initial_conditions,
      rtol = TOL[i],
      atol = TOL[i],
      max_num_steps = MNS[j]
    )

    fit <- rstan::sampling(model, dat, chains = C, refresh = 0)
    runtimes[i, j, ] <- rowSums(get_elapsed_time(fit))
  }
}

# Plotting runtimes for all chains
DF <- NULL
for (i in 1:NTOL) {
  tim <- runtimes[i, , ]
  tim <- as.vector(t(tim))
  tol <- rep(TOL[i], length(tim))
  mns <- rep(MNS, each = C)
  df <- data.frame(tol, mns, tim)
  DF <- rbind(DF, df)
}
DF$tol <- as.factor(DF$tol)
colnames(DF) <- c("tol", "max_num_steps", "runtime")

aes <- aes_string(
  x = "max_num_steps", y = "runtime",
  group = "tol", color = "tol", shape = "tol"
)
plt1 <- ggplot(DF, aes) +
  geom_point() +
  scale_x_log10()

# Plotting runtimes averged over chains
rt <- apply(runtimes, c(1, 2), mean)
tim <- as.vector(t(rt))
mns <- rep(MNS, NTOL)
tol <- rep(TOL, each = L)
df <- data.frame(as.factor(tol), mns, tim)
colnames(df) <- c("tol", "max_num_steps", "runtime")

aes <- aes_string(
  x = "max_num_steps", y = "runtime",
  group = "tol", color = "tol"
)
plt2 <- ggplot(df, aes) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  ylab("average runtime over chains (s)")
