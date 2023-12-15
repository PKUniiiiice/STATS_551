data {
  int<lower=1> n;   // number of data items
  int<lower=1> p;   // number of predictors, contain intercept
  matrix[n, p] X;   // x matrix
}
parameters {
  int<lower=1> K;
  vector[K] p;       
  vector[K] mu[p];
  
  #together Sigma
  vector[p] sigma[K];
  corr_matrix[p] L[K];
  
  
}
transformed parameters {
  real sigma = exp(log_sigma);
}
model {
  y ~ normal(x * beta, sigma);  // likelihood
}
