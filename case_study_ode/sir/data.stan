  int<lower=1> N; // number of time points
  real t[N]; // time points
  int<lower=1> pop_size; // population size
  real<lower=1> I0; // initial number of infected
  real<lower=0> RTOL; // ODE solver relative tolerance
  real<lower=0> ATOL; // ODE solver absolute tolerance
  int<lower=1> MAX_NUM_STEPS; // ODE solver maximum number of steps
