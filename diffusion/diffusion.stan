functions {
#include tridiag.stan
#include be.stan
}

data {
  int N;
  vector[N] x;
  vector[N] y;
  real T;
  real dt;
  real sigma;
}

transformed data {
  real dx = x[2] - x[1];
}

parameters {
  real<lower = 0.0> K; // diffusion constant
}

transformed parameters {
  vector[N] u;
  real t = 0.0;
  
  // Initialize u(t=0,x)
  for(n in 1:N) {
    if(x[n] > 0.5) {
      u[n] = 1.0;
    } else {
      u[n] = 0.0;
    }
  }
  
  // Solve u(t=T,x) using backward Euler method
  u = stan_be(u, dt, dx, T, K);
}

model {
  K ~ normal(0, 1);
  y ~ normal(u, sigma);
}
