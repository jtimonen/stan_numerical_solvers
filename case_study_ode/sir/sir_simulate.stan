functions {
#include functions.stan
}
data {
#include data.stan
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
generated quantities {
  vector[3] x[N] = ode_rk45_tol(SIR, x0, t0, t, RTOL, ATOL, 
    MAX_NUM_STEPS, a0, to_vector({beta, gamma}));
  int y[3, N];
  for(n in 1:N) {
    for(j in 1:3) {
      y[j,n] = neg_binomial_2_rng(x[n][j] + 10*ATOL, phi); 
    }
  }
}
