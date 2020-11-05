data {
  int N;
  int y[N];
}

parameters {
  real log_lambda;
}

transformed parameters {
  real lambda = exp(log_lambda);
}

model {
  y ~ poisson(lambda);
}
