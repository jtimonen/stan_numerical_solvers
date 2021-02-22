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
}

parameters {
  vector<lower=0>[2] theta;
  real<lower=0> sigma;
}

model {
  theta ~ normal(1, 0.3);
  sigma ~ normal(0, 2.0);
  vector[2] y_hat[N] = ode_rk45(derivative_fun, y0, t0, t_eval, a0, theta);
  for(n in 1:N){
    target += normal_lpdf(y_data[n] | y_hat[n], sigma); 
  }
}
