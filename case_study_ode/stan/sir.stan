functions {
#include functions-SIR.stan
}

data {
  int<lower=1> N; // number of observations
  real t[N]; // observation time points
  int y[N]; // measurements of infected
  int<lower=1> P; // number of points where to simulate a solution
  real t_sim[P]; // simulation time points
  int<lower=1> pop_size; // population size
  real<lower=1> I0; // initial number of infected
  real<lower=0> RTOL;
  real<lower=0> ATOL;
  int<lower=1> MAX_NUM_STEPS;
}

transformed data {
  real t0 = 0.0;
  int a0[1] = {pop_size};
  vector[3] x0 = to_vector({pop_size - I0, I0, 0.0}); // {S, I, R}
}

parameters {
  real<lower=0.05, upper = 3> beta;
  real<lower=0.05, upper = 3> gamma;
  real<lower=0.1, upper = 10> phi;
}

transformed parameters {
  vector[2] theta = to_vector({beta, gamma});
}

model {
  vector[3] x[N] = ode_rk45_tol(SIR, x0, t0, t, RTOL, ATOL, 
    MAX_NUM_STEPS, a0, theta);
  // Add small positive number to solution to avoid negative numbers
  for(n in 1:N) {
    target += neg_binomial_2_lpmf(y[n] | x[n] + 10*ATOL, phi);
  }
}

generated quantities {
  vector[3] x_sim[N] = ode_rk45_tol(SIR, x0, t0, t_sim, RTOL, ATOL, 
    MAX_NUM_STEPS, a0, theta);
}
