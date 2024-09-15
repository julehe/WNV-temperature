functions {
  // a linear function
  real f(real x, real alpha, real beta) {
    real out;
    out = exp(alpha-beta*x);
    return out;
  }
}
data {
  int<lower=0> Nobs; // number of observations
  int<lower=0> Nexp; // number of experiments
  int<lower=0> Nindex; // number of experiment and temperature setting combinations
  
  vector[Nobs] dpi; // predictor (days post infection)
  array[Nobs] int<lower=1, upper=Nindex> id; // id of experiment and temperature setting combination
  vector[Nobs] y; // response (trait)
  

  int<lower=0> N_new; // number of temp. points for generated quantities
  vector[N_new] x_new; // temp. points for generated quantities
  int<lower=0> N_dpi_new; // number of days post infection for generated quantities
  vector[N_dpi_new] dpi_new; // days post infection for generated quantities
  vector[Nindex] x_index; // temperature setting of each experiment and temperature setting combination   
  array[Nindex] int<lower=1, upper=Nexp> exp_index; // experiment id of each experiment and temperature setting combination
}
parameters {
  // population-level means
  real mu_alpha;
  real mu_beta;
  
  // between-experiment standard deviations
  vector[Nexp] exp_alpha;
  vector[Nexp] exp_beta;
  
  real<lower=0> sigma_EIP; // EIP standard deviation
  real<lower=0> sigma; // measurement error
}
transformed parameters{
  // experiment-level parameters
  vector[Nexp] alpha;
  vector[Nexp] beta;
  for (n in 1:Nexp) {
    alpha[n] = mu_alpha + exp_alpha[n];
    beta[n] = mu_beta + exp_beta[n];
  }
  
  // likelihood mean for each observation
  vector[Nindex] mu;
  for (n in 1:Nindex) {
    mu[n] = f(x_index[n], alpha[exp_index[n]], beta[exp_index[n]]);
  }
  vector[Nobs] p;
  for (n in 1:Nobs) {
    p[n] = normal_cdf(dpi[n],mu[id[n]],sigma_EIP);
  }
}
model {
  // hyperpriors
  mu_alpha ~ normal(0,10);
  mu_beta ~ normal(0,1);
  
  // priors
  exp_alpha ~ normal(0,0.1);
  exp_beta ~ normal(0,0.01);
  
  sigma_EIP ~ normal(0,50);
  sigma ~ normal(0,1);
  
  // likelihood
  y ~ normal(p, sigma);
}
generated quantities {
  // posterior expected temperature response per experiment
  matrix[Nindex, N_dpi_new] p_new;
  for (i in 1:Nindex){
    for (n in 1:N_dpi_new){
      real temporary = f(x_index[i], alpha[exp_index[i]], beta[exp_index[i]]);
      p_new[i,n] = normal_cdf(dpi_new[n], temporary, sigma_EIP);
    }
  }
  matrix[Nexp, N_new] f_new;
  for (m in exp_index) {
    for (n in 1:N_new) {
      f_new[m, n] = f(x_new[n], alpha[m], beta[m]);
    }
  }
  
  // posterior expected temperature response new experiment
  real alpha_new;
  real beta_new;
  alpha_new = normal_rng(mu_alpha,0.1);
  beta_new = normal_rng(mu_beta,0.01);
  vector[N_new] f_new_spec;
  for (n in 1:N_new) {
    f_new_spec[n] = f(x_new[n], alpha_new, beta_new);
  }
}
