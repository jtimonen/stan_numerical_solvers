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
    
    return {dS_dt, dI_dt, dR_dt};
  }
  
  // Solve the SIR system
  real[,] stan_solve_sir(real[] y0, real[] ts, real[] theta,
      data real[] x_r, data int N, data real rtol, data real atol, 
      data int max_steps) {
    int n_days = num_elements(ts);
    int x_i[1] = { N };
    real f[n_days, 3] = integrate_ode_bdf(stan_sir, y0, 0.0, ts, theta, 
      x_r, x_i, rtol, atol, max_steps);
    return(f);
  }
  
}

data {
  int<lower=1> n_days;  // Number of observations
  real y0[3];           // Initial state (day 0)
  real ts[n_days];      // Observation times
  int N;                // Population size
  int cases[n_days];    // Counts of infected people
  
  // Control parameters
  real<lower=0.0> rtol;
  real<lower=0.0> atol;
  int<lower=100> max_steps;
}

transformed data {
  real x_r[0];
}

parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi;
}

transformed parameters{
  real y_hat[n_days, 3];
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;
    y_hat = stan_solve_sir(y0, ts, theta, x_r, N, rtol, atol, max_steps);
  }
}

model {
  gamma ~ normal(0, 1);
  beta ~ normal(0, 1);
  phi ~ lognormal(1, 1);
  cases ~ neg_binomial_2(y_hat[,2], phi);
}
