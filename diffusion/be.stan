// Backward Euler method for solving the 1D diffusion problem
// 
vector stan_be(vector u_init, real dt, real dx, real T_max, real K){
  real t = 0.0;
  int Nx = num_elements(u_init);
  vector[Nx] u = u_init;
  real K_star = K * dt / (dx^2);
  
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
