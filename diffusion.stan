data {
  int N;
  vector[N] x;
  real L;
  real Tf;
  real dt;
  real sigma;
  vector[N] y;
}

transformed data {
  real dx = x[2] - x[1];
}

parameters {
  real<lower = 0.0> D;
}

transformed parameters {
  vector[N] u;
  real t = 0.0;
  
  for(n in 1:N) {
    if(x[n] > L / 2.0) {
      u[n] = 1.0;
    } else {
      u[n] = 0.0;
    }
  }

  while(t < Tf) {
    vector[N] u_new = u;
  
    for(n in 2:(N - 1)) {
      u_new[n] = dt * (u[n + 1] - 2 * u[n] + u[n - 1]) / (2 * dx) + u[n];
    }
  
    u = u_new;
    t = t + dt;
  }
}

model {
  D ~ normal(0, 1);
  y ~ normal(u, sigma);
}
