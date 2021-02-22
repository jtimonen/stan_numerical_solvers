functions {
#include lv_system.stan
}

data {
  int<lower=1> N;
  real t_eval[N]; // must be increasing
  vector[2] y_data[N];
  vector[2] y0;
  real t0;
}

transformed data {
  int a0[0];
  vector[2*N] y_data_vec;
  for (i in 1:N) {
    y_data_vec[i] = y_data[i][1];
    y_data_vec[N + i] = y_data[i][2];
  }
}

parameters {
  vector<lower=0>[2] theta;
  real<lower=0> sigma;
}

model {
  vector[2*N] y_hat_vec;
  theta ~ normal(1, 0.3);
  sigma ~ normal(0, 2.0);
  vector[2] y_hat[N] = ode_rk45(derivative_fun, y0, t0, t_eval, a0, theta);
  for (i in 1:N) {
    y_hat_vec[i] = y_hat[i][1];
    y_hat_vec[N + i] = y_hat[i][2];
  }
  target += normal_lpdf(y_data_vec | y_hat_vec, sigma);
}
