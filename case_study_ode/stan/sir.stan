functions {
  // SIR system right-hand side
  vector derivative_fun(real t, vector y, data int[] a0, vector theta) {
    int pop_size = a0[1];
    vector[3] dy_dt;
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real beta = theta[1];
    real gamma = theta[2];
    real infection_rate = beta * I * S / pop_size;
    real recovery_rate = gamma * I;
    dy_dt[1] = - infection_rate;
    dy_dt[2] = infection_rate - recovery_rate;
    dy_dt[3] = recovery_rate;
    return dy_dt;
  }
}

data {
  int<lower=1> N; // number of observations
  real t_data[N]; // data time points, must be increasing
  int<lower=1> pop_size; // population size
  real<lower=1> I0; // initial number of infected
  real<lower=0> RTOL;
  real<lower=0> ATOL;
  int<lower=1> MAX_NUM_STEPS;
  
  int y_data[N]; // measurements of infected
}

transformed data {
  real t0 = 0.0;
  int a0[1] = {pop_size};
  vector[3] y0 = to_vector({pop_size - I0, I0, 0.0}); // {S, I, R}
}

parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> phi;
}

transformed parameters {
  real phi_inv = 1.0 / phi;
}

model {
  beta ~ gamma(5, 5);
  gamma ~ gamma(5, 5);
  phi_inv ~ gamma(10, 15);
  vector[3] y_hat[N];
  {
    vector[2] theta = to_vector({beta, gamma});
    y_hat = ode_bdf_tol(derivative_fun, y0, t0, t_data, 
      RTOL, ATOL, MAX_NUM_STEPS, a0, theta);
  }
  
  // Add small positive number to solution to avoid negative numbers
  for(n in 1:N) {
    target += neg_binomial_2_lpmf(y_data[n] | y_hat[n] + 2.0 * ATOL, phi);
  }
}
