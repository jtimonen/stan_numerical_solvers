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
  real sigma = 0.3;
  real t_grid[G];
  t_grid[1] = t0;
  for (j in 1:num_steps) {
    t_grid[j+1] = t0 + j*h;
  }
}

parameters {
  vector<lower=0>[2] theta;
}

transformed parameters {
  vector[2] y_grid_rk4[G] = odeint_rk4(t0, y0, h, num_steps, a0, theta);
  vector[2] y_rk4[N] = interp_1d_cubic(y_grid_rk4, t_grid, t_eval, a0, theta);
  //vector[2] y_ref[N] = ode_rk45(derivative_fun, y0, t0, t_eval, a0, theta);
}

model {
  theta ~ normal(1, 0.3);
  for (n in 1:N) {
    y_data[n] ~ normal(y_rk4[n], sigma);
  }
}
