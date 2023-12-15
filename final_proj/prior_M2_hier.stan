data {
  int<lower=1> N;   // number of data items
  int<lower=1> K;   // number of predictors, contain intercept
  matrix[N, K] x;   // design matrix
  vector[N] y;      // outcome vector
  
  
}
parameters {
  //indicator of predictor
  //int<lower=0> z[K-1];
  
  //parameters use conjugate
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma2;  //variance
  
  //hyperparameters use uniform
  //for beta m0, C0
  vector[K] mu_beta;
  corr_matrix[K] C0;
  
  //for sigma
  real<lower=0> v0;
  real<lower=0> s0;
}
transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
  //vector[K] beta_ind = beta .* append_row(1, z);
}
model {
  //prior of z
  //bernoulli dist
  //z ~ bernoulli(0.5);
  
 // hyperprior
 C0 ~ lkj_corr(1.0);
 
 v0 ~ cauchy(0,1);
 s0 ~ cauchy(0,1);

 //priot
 sigma2 ~ inv_gamma(v0, s0);
 beta ~ multi_normal(mu_beta, sigma2 * C0);

 // likelihood
 y ~ normal(x * beta, sigma);
}
