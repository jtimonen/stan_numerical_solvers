functions {
#include functions-SIR.stan
}

data {
  int<lower=1> N; // number of points where to simulate a solution
  real t[N]; // simulation time points
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
  real<lower=0> beta;
  real<lower=0> gamma;
}

generated quantities {
  
  vector[3] x[N] = ode_rk45_tol(SIR, x0, t0, t, RTOL, ATOL, 
    MAX_NUM_STEPS, a0, to_vector({beta, gamma}));
  
}
