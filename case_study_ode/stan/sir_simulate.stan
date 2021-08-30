functions {
#include functions-SIR.stan
}

data {
  int<lower=1> N; // number of time points
  real t[N]; // time points
  int<lower=1> pop_size; // population size
  real<lower=1> I0; // initial number of infected
  real<lower=0> RTOL;
  real<lower=0> ATOL;
  int<lower=1> MAX_NUM_STEPS;
  
  int<lower=1> S; // number of parameter sets
  vector<lower=0>[2] THETA[S]; // first column beta, second gamma
}

transformed data {
  real t0 = 0.0;
  int a0[1] = {pop_size};
  vector[3] x0 = to_vector({pop_size - I0, I0, 0.0}); // {S, I, R}
}

parameters {
  real dummy;
}

generated quantities {
  vector[3] x[S, N];
  # Solve ODE at t_eval given each parameter set
  for(s in 1:S) {
    x[s, :] = ode_rk45_tol(SIR, x0, t0, t,
    RTOL, ATOL, MAX_NUM_STEPS, a0, THETA[s]);
  }
}
