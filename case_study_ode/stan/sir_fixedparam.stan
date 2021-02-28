functions {
  // SIR system right-hand side
  vector derivative_fun(real t, vector y, data int[] a0, vector theta) {
    int pop_size = a0[1];
    vector[3] dy_dt;
    dy_dt[1] = -theta[1] * y[2] * y[1] / pop_size;
    dy_dt[2] =  theta[1] * y[2] * y[3] / pop_size - theta[2] * y[2];
    dy_dt[3] =  theta[2] * y[2];
    return dy_dt;
  }
}

data {
  int<lower=1> N; // number of observations
  real t_eval[N]; // output time points, must be increasing
  int<lower=1> pop_size; // population size
  real<lower=1> I0; // initial number of infected
  real<lower=0> RTOL;
  real<lower=0> ATOL;
  int<lower=1> MAX_NUM_STEPS;
  
  int<lower=1> S; // number of parameter sets
  vector<lower=0>[2] THETA[S]; // first column beta, second gamma
}

transformed data {
  real t0 = 0.0;
  int a0[1] = {pop_size};
  vector[3] y0 = to_vector({pop_size - I0, I0, 0.0}); // {S, I, R}
}

parameters {
  real dummy; //ignore this
}

generated quantities {
  vector[2] y_hat[S, N];
  for(s in 1:S) {
    y_hat[s, :] = ode_rk45_tol(derivative_fun, y0, t0, t_eval, 
      RTOL, ATOL, MAX_NUM_STEPS, a0, THETA[s]);
  }
}
