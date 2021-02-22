// Get indices of left points
int[] indices_left(data real[] x, data real[] x_out) {
  int idx = 1;
  int N_out = size(x_out);
  int inds[N_out];
  for (j in 1:N_out) {
    while(x[idx+1] < x_out[j]) idx += 1;
    inds[j] = idx;
  }
  return(inds);
}

// Get indices of right points
int[] indices_right(data real[] x, data real[] x_out) {
  int idx = 1;
  int N_out = size(x_out);
  int inds[N_out];
  for (j in 1:N_out) {
    while(x[idx] < x_out[j]) idx += 1;
    inds[j] = idx;
  }
  return(inds);
}

// Get breaks
int[,] interval_indices(data real[] x_grid, data real[] x_eval) {
  int idx = 1;
  int idx_prev = 0;
  int N_grid = size(x_grid);
  int N_eval = size(x_eval);
  int inds[N_grid-1, 2];
  for (j in 1:(N_grid-1)) {
    while(x_eval[idx] < x_grid[j+1] && idx < N_eval) {
      idx += 1;
    }
    inds[j,] = {idx_prev+1, idx};
    idx_prev = idx;
  }
  return(inds);
}

// Cubic interpolation using Hermite splines
vector[] interp_1d_cubic(vector[] y, data real[] x, data real[] x_out,
    int[,] interval_idx, int[] a0, vector theta){
  int N_in = size(x);
  int N_out = size(x_out);
  int D = num_elements(y[1]);
  vector[N_out] x_out_vec = to_vector(x_out);
  vector[D] f_left;
  vector[D] f_right;
  vector[N_out] y_out[D];
  for (i in 1:(N_in-1)) {
    
    int i0 = interval_idx[i,1];
    int i1 = interval_idx[i,2];
    int L = i1 - i0 + 1;
    
    if (L > 0) {
    
      // Hermite basis functions
      real dx = x[i+1] - x[i];
      vector[L] w = (x_out_vec[i0:i1] - x[i]) / dx;
      vector[L] w2 = w .* w;
      vector[L] w3 = w2 .* w;
      vector[L] h00 = 2.0 * w3 - 3.0 * w2 + 1.0;
      vector[L] h10 = w3 - 2.0 * w2 + w;
      vector[L] h01 = -2.0 * w3 + 3.0 * w2;
      vector[L] h11 = w3 - w2;
      
      // Evaluate derivatives
      f_left = derivative_fun(0.0, y[i], a0, theta);
      f_right = derivative_fun(0.0, y[i+1], a0, theta);
    
      // Compute interpolation
      for(d in 1:D) {
        y_out[d][i0:i1] = h00*y[i,d] + dx*h10*f_left[d] + 
        h01*y[i+1, d] + dx*h11*f_right[d];
      }
    }
  }
  return(y_out);
}

// RK4 method
vector[] odeint_rk4(real t0, vector y0, data real h, data int num_steps,
    int[] a0, vector theta){

  int d = num_elements(y0);
  vector[d] y[num_steps+1];
  vector[d] k1; vector[d] k2; vector[d] k3; vector[d] k4;
  real t = t0;
  y[1] = y0;
  
  // Integrate at grid of time points
  for(i in 1:num_steps){
    k1 = h * derivative_fun(t          , y[i]           , a0, theta);
    k2 = h * derivative_fun(t + 0.5 * h, y[i] + 0.5 * k1, a0, theta);
    k3 = h * derivative_fun(t + 0.5 * h, y[i] + 0.5 * k2, a0, theta);
    k4 = h * derivative_fun(t + h      , y[i] + k3      , a0, theta);
    y[i+1] = y[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    t = t + h;
  }
  
  return(y);
}

