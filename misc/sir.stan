functions {
  // SIR system right-hand side
  real[] stan_sir(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real N = x_i[1];
    
    real beta = theta[1];
    real gamma = theta[2];
    
    real dS_dt = -beta * I * S / N;
    real dI_dt =  beta * I * S / N - gamma * I;
    real dR_dt =  gamma * I;
    
    return {dS_dt, dI_dt, dR_dt};
  }
  
}

data {
  
}
