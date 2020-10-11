functions {
  // Solve a symmetric tridiagonal linear system Ax = d with constant
  // secondary diagonals
  // 
  // a = the constant value of both secondary diagonals of matrix A
  // b = diagonal of matrix A
  vector stan_solve_tridiag_be(real a, vector b, vector d){
    int n = num_elements(b);
    vector[n] x = rep_vector(0.0, n);
    real w;
    int idx;
    vector[n] bb = b;
    vector[n] dd = d;
    
    // Forward sweep
    for (i in 2:n) {
      w = a / bb[i - 1];
      bb[i] = bb[i] - w*a;
      dd[i] = dd[i] - w*dd[i - 1];
    }
    
    // Back substitution
    x[n] = dd[n]/bb[n];
    idx = n - 1;
    while(idx > 0) {
      x[idx] = (dd[idx] - a*x[idx + 1]) / bb[idx];
      idx = idx - 1;
    }
    return(x);
  }
  
  // Backward Euler method for solving the 1D diffusion problem
  //
  // u_init, initial conditions
  // dt, timestep
  // dx, spatial discretization
  // T_max, max time
  // K, diffusion constant
  // ul, left boundary condition
  // ur, right boundary condition
  vector stan_be(vector u_init, real dt, real dx, real T_max, real K, real ul, real ur){
    int Nx = num_elements(u_init);
    vector[Nx] u = u_init;
    real K_star = K * dt / (dx^2); // this value defines some properties
    real t = 0.0;
    
    // Create the diagonal of the tridiagonal matrix A
    vector[Nx] A_diag = rep_vector(1.0 + 2.0 * K_star, Nx);
  
    // Iterate time step
    while(t < T_max) {
      vector[Nx] b = u;
      b[1] += ul * K_star;
      b[Nx] += ur * K_star;
      // Update u and t
      u = stan_solve_tridiag_be(-K_star, A_diag, b);
      t = t + dt;
    }
    return(u);
  }
}

data {
  int N;
  vector[N] x;
  vector[N] y;
  real T;
  real dt;
  real sigma;
  real ul;
  real ur;
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
  u = stan_be(u, dt, dx, T, K, ul, ur);
}

model {
  K ~ normal(0, 1.0);
  y ~ normal(u, sigma);
}
