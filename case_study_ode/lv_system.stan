  // Lotka-Volterra system
  vector derivative_fun(real t, vector y, data int[] a0, vector theta) {
    vector[2] dydt;
    dydt[1] = theta[1]*y[1] - y[1]*y[2];
    dydt[2] = y[1]*y[2] - theta[2]*y[2];
    return dydt;
  }
