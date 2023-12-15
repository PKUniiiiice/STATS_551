data {
  int<lower=1> N;   // number of data items
  int<lower=1> K;   // number of predictors, contain intercept
  matrix[N, K] x;   // design matrix
  vector[N] y;      // outcome vector
  
  vector[K] b0;
}
parameters {
  real<lower=0> g;
  
  //paras
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // standard deviation
}
model {
 // prior
 g ~ inv_gamma(1.0/2.0, N/2.0);
 beta ~ multi_normal(b0, square(sigma) * g * inverse_spd(x' * x));
 target += -2 * log(sigma);
 
 // likelihood
 y ~ normal(x * beta, sigma);
}
