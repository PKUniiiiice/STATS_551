data {
  int<lower=1> N;   // number of data items
  int<lower=1> K;   // number of predictors, contain intercept
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
}
parameters {
  vector[K] beta;       // coefficients for predictors
  real log_sigma;  // Jeffreys' prior on log is uniform
}
transformed parameters {
  real sigma = exp(log_sigma);
}
model {
  y ~ normal(x * beta, sigma);  // likelihood
}
