functions {
  real approx_exp(real dx, real x) {
    real y = 1.0;
    real tx = 0.0;

    while(fabs(tx + dx) < abs(x)) {
      y += dx * y;
      tx += dx;
    }
    y += (x - tx) * y;
    
    return y;
  }
}

data {
  int N;
  int y[N];
  real dx;
}

parameters {
  real log_lambda;
}

transformed parameters {
  real lambda = approx_exp(dx, log_lambda);
}

model {
  y ~ poisson(lambda);
}

generated quantities {
  real ref_lambda = exp(log_lambda);
  real log_ratio = poisson_lpmf(y | ref_lambda) - poisson_lpmf(y | lambda);
}
