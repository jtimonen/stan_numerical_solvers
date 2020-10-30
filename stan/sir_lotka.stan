functions {
  // SIR system right-hand side
  row_vector stan_sir(real t, row_vector y,
                      real alpha, real beta,
                      real gamma, real delta) {
    real u = y[1];
    real v = y[2];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;

    return [ du_dt, dv_dt ];
  }
  
  // Solve the SIR system
  matrix stan_solve_sir(row_vector initial_conditions, vector ts,
                        real alpha, real beta,
                        real gamma, real delta,
                        data real dt) {
    int N = num_elements(ts);
    
    vector[N] I;
    
    matrix[N, 2] mu;
    real t = ts[1];
    // We're not saving the initial conditions -- only the output times
    row_vector[2] y1 = initial_conditions;
    for(i in 1:N) {
      // Break from loop when we have integrated past an output
      while(t < ts[i]) {
        // Integrate with explicit midpoint method
        y1 = y1 + dt * stan_sir(t, y1, alpha, beta, gamma, delta);
        // Integrate with explicit midpoint method
        //row_vector[2] yh = y1 + (dt / 2.0) * stan_sir(t, y1, alpha, beta, gamma, delta);
        //y1 = y1 + dt * stan_sir(t + dt / 2.0, yh, alpha, beta, gamma, delta);
        t += dt;
      }
      // Use euler method to figure out what time point was
      mu[i, ] = y1 - (t - ts[i]) * stan_sir(t, y1, alpha, beta, gamma, delta);
    }
    
    return mu;
  }
}
  
data {
  int<lower=1> N;             // Number of observations
  vector[N] ts;                 // Observation times
  real<lower = 0.0> y[N, 2];                   // Counts of infected people
  //row_vector[2] initial_conditions; // Initial state (day 0)
  
  real<lower=0.0> dt;         // Timestep
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> delta;
  row_vector<lower = 0.0>[2] initial_conditions;
  real<lower = 0.0> sigma[2];
}

transformed parameters {
  matrix[N, 2] mu = stan_solve_sir(initial_conditions, ts,
                               alpha, beta, gamma, delta, dt);
}

model {
  alpha ~ normal(1, 0.5);
  gamma ~ normal(1, 0.5);
  beta ~ normal(0.05, 0.05);
  delta ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  initial_conditions ~ lognormal(log(10), 1);
  for (k in 1:2) {
    y[ , k] ~ lognormal(log(mu[, k]), sigma[k]);
  }
}
