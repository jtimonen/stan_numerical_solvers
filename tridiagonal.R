N = 4

b0 = rnorm(N)
A = -2 * diag(N)

for (i in 1:(N - 1)) {
  A[i, i + 1] = 1.0;
  A[i + 1, i] = 1.0;
}

dt = 0.23
dx = 0.12
D = 0.25

alpha = D * dt / dx^2

# Do the solve with a dense matrix
r1 = solve(diag(N) - alpha * A, b0)

# 'rl' is my placeholder name for the offdiagonal elements of the matrix alpha * A
rl = -alpha

b = b0
d = 1 - diag(alpha * A)
for (i in 2:N) {
  b[i] = b[i] - (rl / d[i - 1]) * b[i - 1]
  d[i] = d[i] - rl^2 / d[i - 1]
}

for (i in (N - 1):1) {
  b[i] = b[i] - (rl / d[i + 1]) * b[i + 1]
}

# Here's the solution the fast way!
r2 = b / d
