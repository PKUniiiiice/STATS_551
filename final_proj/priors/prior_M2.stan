data {
  int<lower=1> N;   // number of data items
  int<lower=1> K;   // number of predictors, contain intercept
  matrix[N, K] x;   // design matrix
  vector[N] y;      // outcome vector
  
  // parameters for the priors
  // for beta
  vector[K] m0;
  matrix[K, K] C0;
  //for sigma^2
  real<lower=0> v0;
  real<lower=0> s0;
  
}
parameters {
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma2;  //variance
}
transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
}
model {
 // prior
 sigma2 ~ inv_gamma(v0, s0);
 beta ~ multi_normal(m0, sigma2 * C0);

 // likelihood
 y ~ normal(x * beta, sigma);
}
