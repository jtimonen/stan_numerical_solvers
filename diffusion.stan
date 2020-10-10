functions {
  // Solve a tridiagonal linear system Ax = d
  // 
  // a = lower diagonal of A
  // b = diagonal of matrix A
  // c = upper diagonal of A
  vector stan_solve_tridiag(vector a, vector b, vector c, vector d){
    int n = num_elements(b);
    vector[n] x = rep_vector(0.0, n);
    real w;
    int idx;
    vector[n-1] aa = a;
    vector[n] bb = b;
    vector[n-1] cc = c;
    vector[n] dd = d;
  
    // Forward sweep
    for (i in 2:n) {
      w = aa[i - 1] / bb[i - 1];
      bb[i] = bb[i] - w*cc[i - 1];
      dd[i] = dd[i] - w*dd[i - 1];
    }
  
    // Back substitution
    x[n] = dd[n]/bb[n];
    idx = n - 1;
    while(idx > 0) {
      x[idx] = (dd[idx] - cc[idx]*x[idx + 1]) / bb[idx];
      idx = idx - 1;
    }
    return(x);
  }
}

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
