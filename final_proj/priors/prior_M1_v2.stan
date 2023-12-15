data {
  int<lower=1> N;   // number of data items
  int<lower=1> K;   // number of predictors, contain intercept
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
}
parameters {
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // standard deviation
}
model {
  //improper prior
  target += -2 * log(sigma); 
  // likelihood
  y ~ normal(x * beta, sigma);  
}
