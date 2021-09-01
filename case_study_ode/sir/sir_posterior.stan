functions {
#include functions.stan
}
data {
#include data.stan
  int y[N]; // measurements of infected
}
transformed data {
#include tdata.stan
}
parameters {
#include params.stan
}
transformed parameters {
#include tparams.stan
}
model {
#include prior.stan
  vector[3] x[N] = ode_rk45_tol(SIR, x0, t0, t, RTOL, ATOL, 
        MAX_NUM_STEPS, a0, to_vector({beta, gamma}));
  for(n in 1:N) {
    y[n] ~ neg_binomial_2(x[n] + 10*ATOL, phi);
  }
}
