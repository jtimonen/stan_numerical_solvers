  // SIR system right-hand side
  vector SIR(real t, vector y, data int[] a0, vector theta) {
    vector[3] dx_dt;
    int pop_size = a0[1];
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real beta = theta[1];
    real gamma = theta[2];
    real infection_rate = beta * I * S / pop_size;
    real recovery_rate = gamma * I;
    dx_dt[1] = - infection_rate;
    dx_dt[2] = infection_rate - recovery_rate;
    dx_dt[3] = recovery_rate;
    return dx_dt;
  }
