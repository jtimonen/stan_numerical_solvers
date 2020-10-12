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
  
  // Solve the SIR system
  // TODO: fix
  // TODO: should we use ode_bdf_tol?
  vector[] stan_integrate_ode(y0, t0, ts, theta, int N,
  real rtol, real atol, int max_steps) {
    vector[n_days] f[3] = integrate_ode_bdf(stan_sir, y0, t0, ts, theta, 
      [], { N }, rtol, atol, max_steps);
    return(f);
  }
  
}

data {
  int<lower=1> n_days;
  real y0[3];
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
}

parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_days, 3];
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}


