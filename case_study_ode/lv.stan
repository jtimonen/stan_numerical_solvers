functions {
#include lv_system.stan
#include odeint.stan
}

data {
  int<lower=1> N;
  real t_eval[N]; // must be increasing
  vector[2] y_data[N];
  vector[2] y0;
  real t0;
  int<lower=1> num_steps; // number of steps
  real<lower=0> h; // step size
}

transformed data {
  int a0[0];
  int G = num_steps + 1;
  real t_grid[G];
  t_grid[1] = t0;
  for (j in 1:num_steps) {
    t_grid[j+1] = t0 + j*h;
  }
}

parameters {
  vector<lower=0>[2] theta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[2] y_inf[N];
  {
    vector[2] y_grid[G] = odeint_rk4(t0, y0, h, num_steps, a0, theta);
    y_inf = interp_1d_cubic(y_grid, t_grid, t_eval, a0, theta);
  }
}

model {
  theta ~ normal(1, 0.3);
  sigma ~ inv_gamma(5, 5);
  for (n in 1:N) target += normal_lpdf(y_data[n] | y_inf[n], sigma);
}

generated quantities {
  real log_ratio;
  {
    real log_lik_inf = 0.0;
    real log_lik_ref = 0.0;
    vector[2] y_ref[N] = ode_rk45(derivative_fun, y0, t0, t_eval, a0, theta);
    for (n in 1:N) {
      log_lik_inf += normal_lpdf(y_data[n] | y_inf[n], sigma);
      log_lik_ref += normal_lpdf(y_data[n] | y_ref[n], sigma);
    }
    log_ratio = log_lik_inf - log_lik_ref;
  }
}
