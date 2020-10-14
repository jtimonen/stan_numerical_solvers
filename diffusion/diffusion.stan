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
  vector solve_pde(real dt, int Nx, real K, real T_meas, vector x_meas) {
    real L = 1.0; // length of rod
    real ul = 1.0; // left boundary condition
    real ur = 0.0; // right boundary condition

    real dx = L / (Nx + 1);
    real K_star = K * dt / (dx^2);
    real t = 0.0;

    vector[Nx] u = rep_vector(0, Nx);
    vector[Nx] u_prev = u;

    vector[rows(x_meas)] solution;

    // Create the diagonal of the tridiagonal matrix A
    vector[Nx] A_diag = rep_vector(1.0 + 2.0 * K_star, Nx);

    // Iterate time step
    while(t < T_meas) {
      vector[Nx] b = u;
      u_prev = u;
  
      b[1] += ul * K_star;
      b[Nx] += ur * K_star;
      // Update u and t
      u = stan_solve_tridiag_be(-K_star, A_diag, b);
      t = t + dt;
    }

    // Use linear interpolation to get solution at T_max not a multiple of dt
    if(T_meas < t) {
      real alpha = (t - T_meas) / dt;
      u = alpha * u_prev + (1.0 - alpha) * u;
    }
    
    // Use linear interpolation to get solution at measurement points
    {
      int i = 1;
      int j = 0;
      while(i <= rows(x_meas)) {
        if(x_meas[i] < 0.0) {
          solution[i] = ul;
          i += 1;
        } else if(x_meas[i] >= L) {
          solution[i] = ur;
          i += 1;
        } else if(j == 0) {
          if(x_meas[i] < dx) {
            real alpha = (dx - x_meas[i]) / dx;
            solution[i] = alpha * ul + (1.0 - alpha) * u[1];
            i += 1;
          } else {
            j += 1;
          }
        } else if(j + 1 == Nx + 1) {
          if(x_meas[i] < L) {
            real alpha = (L - x_meas[i]) / dx;
            solution[i] = alpha * u[Nx] + (1.0 - alpha) * ur;
            i += 1;
          } else {
            j += 1;
          }
        } else {
          if(x_meas[i] >= j * dx && x_meas[i] < (j + 1) * dx) {
            real alpha = ((j + 1) * dx - x_meas[i]) / dx;
            solution[i] = alpha * u[j] + (1.0 - alpha) * u[j + 1];
            i += 1;
          } else {
            j += 1;
          }
        }
      }
    }

    return(solution);
  }
}

data {
  real dt;
  int Nx;
  int N_meas;
  real T_meas;
  vector[N_meas] x_meas;
  
  vector[N_meas] y;
}

parameters {
  real<lower = 0.0> K;     // diffusion constant
  real<lower = 0.0> sigma; // noise magnitude
}

transformed parameters {
  vector[N_meas] mu = solve_pde(dt, Nx, K, T_meas, x_meas);
}

model {
  sigma ~ normal(0, 1.0);
  K ~ normal(0, 1.0);
  y ~ normal(mu, sigma);
}
