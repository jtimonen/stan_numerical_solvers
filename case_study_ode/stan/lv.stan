functions {
#include lv_system.stan
#include odeint.stan
}

data {
  int<lower=1> N;
  real t_eval[N]; // must be increasing
  vector[N] y_data[2];
  vector[2] y0;
  real t0;
  int<lower=1> num_steps; // number of steps
  real<lower=0> h; // step size
  int<lower=0, upper=1> SOLVER;
  real<lower=0> ATOL;
  real<lower=0> RTOL;
  int<lower=10> max_num_steps;
}

transformed data {
  int a0[0];
  int G = num_steps + 1;
  int interval_idx[N-1, 2];
  real t_grid[G];
  {
    t_grid[1] = t0;
    for (j in 1:num_steps) {
      t_grid[j+1] = t0 + j*h;
    }
    interval_idx = interval_indices(t_grid, t_eval);
  }
}

parameters {
  vector<lower=0>[2] theta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] y_hat[2];
  if (SOLVER == 0) {
    vector[2] y_grid[G] = odeint_rk4(t0, y0, h, num_steps, a0, theta);
    y_hat = interp_1d_cubic(y_grid, t_grid, t_eval, interval_idx, a0, theta);
  } else {
    //y_hat = ode_rk45_tol(derivative_fun, y0, t0, t_eval, 
    //  RTOL, ATOL, max_num_steps, a0, theta);
  }
}

model {
  theta ~ normal(1, 0.3);
  sigma ~ inv_gamma(5, 5);
  for(j in 1:2){
    target += normal_lpdf(y_data[j] | y_hat[j], sigma); 
  }
}
