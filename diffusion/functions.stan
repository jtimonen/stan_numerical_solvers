// Solve a tridiagonal linear system Ax = d
// 
// a = lower diagonal of A
// b = diagonal of matrix A
// c = upper diagonal of A
vector stan_solve_tridiag(vector a, vector b, vector c, vector d){
  int n = num_elements(b);
  vector[n] x = rep_vector(0.0, n);
  real w;
  int idx;
  vector[n] bb = b;
  vector[n] dd = d;
  
  // Forward sweep
  for (i in 2:n) {
    w = a[i - 1] / bb[i - 1];
    bb[i] = bb[i] - w*c[i - 1];
    dd[i] = dd[i] - w*dd[i - 1];
  }
  
  // Back substitution
  x[n] = dd[n]/bb[n];
  idx = n - 1;
  while(idx > 0) {
    x[idx] = (dd[idx] - c[idx]*x[idx + 1]) / bb[idx];
    idx = idx - 1;
  }
  return(x);
}

// Solve a symmetric tridiagonal linear system Ax = d
// 
// a = lower/upper diagonal of A
// b = diagonal of matrix A
vector stan_solve_tridiag_sym(vector a, vector b, vector d){
  int n = num_elements(b);
  vector[n] x = rep_vector(0.0, n);
  real w;
  int idx;
  vector[n] bb = b;
  vector[n] dd = d;
  
  // Forward sweep
  for (i in 2:n) {
    w = a[i - 1] / bb[i - 1];
    bb[i] = bb[i] - w*a[i - 1];
    dd[i] = dd[i] - w*dd[i - 1];
  }
  
  // Back substitution
  x[n] = dd[n]/bb[n];
  idx = n - 1;
  while(idx > 0) {
    x[idx] = (dd[idx] - a[idx]*x[idx + 1]) / bb[idx];
    idx = idx - 1;
  }
  return(x);
}

// Solve a symmetric tridiagonal linear system Ax = d with constant band
// 
// a = the constant value of lower/upper diagonal of A
// b = diagonal of matrix A
vector stan_solve_tridiag_be(real a, vector b, vector d){
  int n = num_elements(b);
  vector[n] x = rep_vector(0.0, n);
  real w;
  int idx;
  vector[n] bb = b;
  vector[n] dd = d;
  
  // Forward sweep
  for (i in 2:n) {
    w = a / bb[i - 1];
    bb[i] = bb[i] - w*a;
    dd[i] = dd[i] - w*dd[i - 1];
  }
  
  // Back substitution
  x[n] = dd[n]/bb[n];
  idx = n - 1;
  while(idx > 0) {
    x[idx] = (dd[idx] - a*x[idx + 1]) / bb[idx];
    idx = idx - 1;
  }
  return(x);
}

