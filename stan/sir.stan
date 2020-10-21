functions {
  // SIR system right-hand side
  real[] stan_sir(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real N = x_i[1];
    
    real beta = theta[1];
    real gamma = theta[2];
    
    real dS_dt = -beta * I * S / N;
    real dI_dt =  beta * I * S / N - gamma * I;
    real dR_dt =  gamma * I;
    
    return { dS_dt, dI_dt, dR_dt };
  }
  
  // Solve the SIR system
  vector stan_solve_sir(real[] y0, real[] ts, real[] theta,
      data real[] x_r, data int M, data real rtol, data real atol,
      data real max_num_steps) {
    int n_days = num_elements(ts);
    int x_i[1] = { M };
    real f[n_days, 3] = integrate_ode_rk45(stan_sir, y0, 0.0, ts, theta, 
      x_r, x_i, rtol, atol, max_num_steps);
    return(to_vector(f[, 2]));
  }
  
}

data {
  int<lower=1> N;             // Number of observations
  real ts[N];                 // Observation times
  int y[N];               // Counts of infected people
  int M;                      // Population size
  real initial_conditions[3]; // Initial state (day 0)
  
  // Control parameters
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
  vector[N] mu = stan_solve_sir(initial_conditions, ts, { beta, gamma },
                                 x_r, M, rtol, atol, max_num_steps);
}

model {
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  phi ~ lognormal(1, 1);
  
  // Add small positive number to solution to avoid negative numbers
  y ~ neg_binomial_2(mu + 2.0 * atol, phi);
}
