data {
  int<lower=1> N;   // number of data items
  int<lower=1> K;   // number of predictors, contain intercept
  matrix[N, K] x;   // design matrix
  vector[N] y;      // outcome vector
  
  // parameters for the priors
  // for beta
  vector[K] b0;
  real<lower=0> g;
  
}
parameters {
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // standard deviation
}
transformed parameters {
  real sigma2 = square(sigma);
}
model {
 // prior
 beta ~ normal(b0, sqrt(g)*sigma2);
 target += -2 * log(sigma);
 // likelihood
 y ~ normal(x * beta, sigma);
}
