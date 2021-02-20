
// Cubic interpolation using Hermite splines
vector[] interp_1d_cubic(vector[] y, data real[] x, data real[] x_out,
    int[] a0, vector theta){
  int left = 1;
  int right = 1;
  real h = 0.0;
  real w = 0.0;
  int N_in = size(x);
  int N_out = size(x_out);
  int D = num_elements(y[1]);
  vector[D] f_left;
  vector[D] f_right;
  real h00; real h10; real h01; real h11;
  vector[D] y_out[N_out];
  for (j in 1:N_out) {
    
    // Find left and right point indices
    while(x[right] < x_out[j]) {
      right = right + 1;
    }
    while(x[left+1] < x_out[j]) {
      left = left + 1;
    }
    
    // Evaluate derivatives
    f_left = derivative_fun(0.0, y[left], a0, theta);
    f_right = derivative_fun(0.0, y[right], a0, theta);
    
    // Hermite basis functions
    h = x[right] - x[left];
    w = (x_out[j] - x[left]) / h;
    h00 = 2.0 * w^3 - 3.0 * w^2 + 1.0;
    h10 = w^3 - 2.0 * w^2 + w;
    h01 = -2.0 * w^3 + 3.0 * w^2;
    h11 = w^3 - w^2;
    
    // Compute interpolation
    y_out[j] = h00 * y[left] + h10 * h * f_left + h01 * y[right] +
        h11 * h * f_right;
  }
  return(y_out);
}

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

