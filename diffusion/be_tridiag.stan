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
vector stan_be(vector u_init, real dt, real dx, real T_max, real K){
  int Nx = num_elements(u_init);
  vector[Nx] u = u_init;
  real K_star = K * dt / (dx^2); // this value defines some properties
  real t = 0.0;
  
  // Create the diagonal of the tridiagonal matrix A
  vector[Nx] A_diag = rep_vector(1.0 + 2.0*K_star, Nx);
  A_diag[1] = 1.0 + K_star;
  A_diag[Nx] = 1.0 + K_star;
  
  // Iterate time step
  while(t < T_max) {
    
    // Update u and t
    u = stan_solve_tridiag_be(-K_star, A_diag, u);
    t = t + dt;
    
    // Ensure that boundary values don't change
    u[1] = u_init[1];
    u[Nx] = u_init[Nx];
  }
  return(u);
}
