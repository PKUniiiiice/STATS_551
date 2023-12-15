data {
  //sample size 442
  int<lower=1> n;
  //number of S4.new (group levels)
  int<lower=1> n_gp;
  //number of fix 
  int<lower=1> p_fix;
  //number of rand slope (contain intercept)
  int<lower=1> p_rand;
  
  //response
  real y[n];
  //matrix of fix age sex(as factor) bmi bp(map) S6(glu)
  matrix[n, p_fix] X1;
  //matrix of random slope S1(tc), S2.new.new(tc/ldl--avoid high corr), S3(hdl), S5(ltg)
  matrix[n, p_rand] X2;
  //which group
  int group[n];
  
  //g-prior for fix
  int<lower=0> g;
}
transformed data {
  //sigma of beta_fix
  matrix[p_fix, p_fix] Sigma_fix;
  Sigma_fix = inverse_spd(X1' * X1);
}
parameters {
  //9 in total
  //priors
  vector[p_fix] beta_fix; 
  vector[p_rand] gamma_rand[n_gp];
  real<lower=0> sigma_gp[n_gp];
  
  //hyper priors
  
  //mean of normal beta fix, var is already def in transfored data
  vector[p_fix] mu_beta_fix;

  //mean of normal gamma
  vector[p_rand] mu_gamma[n_gp];
  //LKJ prior for Sigma of gamma
  vector<lower=0>[p_rand] sigma_gamma;
  corr_matrix[p_rand] L;
  
  //mean of sigma_gp
  real mu_sigma_gp;
  real<lower=0> tau_sigma;
}
model {
  //prior
  beta_fix ~ multi_normal(mu_beta_fix, g*Sigma_fix);
  
  gamma_rand ~ multi_normal(mu_gamma, diag_matrix(sigma_gamma) * L * diag_matrix(sigma_gamma));
  
  for (i in 1:n_gp) sigma_gp[i] ~ lognormal(mu_sigma_gp, tau_sigma);
  
  //hyper prior
  //mu_beta_fix unif, mu_gamma unif, mu_sigma_gp unif
  sigma_gamma ~ cauchy(0,1);
  tau_sigma ~ cauchy(0,1);
  L ~ lkj_corr(1.0);
  
  //likelihood
  for (j in 1:n) y[j] ~ normal(X1[j,]*beta_fix + X2[j,]*gamma_rand[group[j]], sigma_gp[group[j]]);
  
}





