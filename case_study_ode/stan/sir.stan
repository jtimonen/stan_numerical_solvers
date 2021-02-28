functions {

  
  // Solve the SIR system
  vector stan_solve_sir(data real[] ts, real[] theta,
      data real[] x_r,
      data real rtol, data real atol, data int max_num_steps) {
    int N = num_elements(ts);
    int M = 1000; // population size
    int I0 = 20; // number of infected on day 0
    int x_i[1] = { M }; // population size
    real y0[3] = { M - I0, I0, 0.0 }; // S, I, R on day 0
    real f[N, 3] = integrate_ode_rk45(stan_sir, y0, 0.0, ts, theta, 
      x_r, x_i, rtol, atol, max_num_steps);
    return(to_vector(f[, 2]));
  }
}

data {
  int<lower=1> N;         // Number of observations
  real t_data[N];         // Observation times
  int y_data[N];          // Counts of infected people
  real<lower=0.0> rtol;
  real<lower=0.0> atol;
  int<lower=1> max_num_steps;
}

transformed data {
  real x_r[0];
}

parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> phi;
}

transformed parameters{
  vector[N] mu = stan_solve_sir(t_data, { beta, gamma },
                                 x_r, rtol, atol, max_num_steps);
}

model {
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  phi ~ lognormal(1, 1);
  
  // Add small positive number to solution to avoid negative numbers
  y_data ~ neg_binomial_2(mu + 2.0 * atol, phi);
}
