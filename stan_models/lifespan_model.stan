functions {
  // a modified linear function
  real linear(real x, real Tmax, real beta) {
    real out;
    real alpha = Tmax*beta;
    if (alpha - beta * x>0) {
      out = alpha - beta * x;
    }else{
      out = 0;
    }
    out = out;
    return out;
  }
}
data {
  int<lower=0> Nspecies; // number of species
  int<lower=0> Nexp; // number of experiments
  int<lower=0> Nobs; // number of observations
  
  vector[Nobs] x; // predictor (temperature)
  vector[Nobs] y; // response (trait performance)
  
  array[Nobs] int<lower=1, upper=Nspecies> species_id; // encodes to which species each observation belongs 
  array[Nobs] int<lower=1, upper=Nexp> exp_id; // encodes to which experiment each observation belongs

  int<lower=0> N_new; // number of temp. points for generated quantities
  vector[N_new] x_new; // temp. points for generated quantities
}
parameters {
  // population-level means
  real mu_beta; 
  real mu_Tmax; 
  
  // between-species standard deviations
  real<lower=0> sigma_beta; 
  real<lower=0> sigma_Tmax; 
  
  // (non-)centered species-level parameters
  vector[Nspecies] Tmax_raw;
  vector<lower=0>[Nspecies] beta;
  
  // between-experiment standard deviations
  real<lower=0> sigma_exp_Tmax;
  real<lower=0> sigma_exp_beta;
  
  // non-centered experiment-level parameters (deviations from species mean)
  vector[Nexp] exp_Tmax;
  vector[Nexp] exp_beta;
  
  real<lower=0> sigma; // measurement error
}
transformed parameters{
  // centered species-level parameters
  vector[Nspecies] Tmax;
  Tmax = Tmax_raw * sigma_Tmax + mu_Tmax;
  
  // species-level intercept 
  vector[Nspecies] alpha;
  for (n in 1:Nspecies) {
    alpha[n] = Tmax[n]*beta[n];
  }
  
  // centered experiment-level parameters (deviations from species mean)
  vector[Nexp] exp_Tmax_sc = exp_Tmax*sigma_exp_Tmax;
  vector[Nexp] exp_beta_sc = exp(exp_beta*sigma_exp_beta);
  
  // likelihood mean for each observation
  vector[Nobs] mu;
  for (n in 1:Nobs) {
    mu[n] = linear(x[n], Tmax[species_id[n]]+exp_Tmax_sc[exp_id[n]], beta[species_id[n]]*exp_beta_sc[exp_id[n]]);
  }
}
model {
  // hyperpriors
  mu_Tmax ~ normal(35,10); 
  mu_beta ~ normal(1.5,1);
  
  sigma_Tmax ~ gamma(4.86,1.88);
  sigma_beta ~ normal(0,1);
  
  sigma_exp_Tmax ~ gamma(18.93,10.07);
  sigma_exp_beta ~ normal(0,1);
  
  // priors
  exp_Tmax ~ std_normal();
  exp_beta ~ std_normal();
  
  beta ~ lognormal(mu_beta,sigma_beta);
  Tmax_raw ~ normal(0,1);
  
  sigma ~ normal(0,10);
  
  // likelihood
  y ~ normal(mu, sigma);
}
generated quantities {
  // posterior predicition replicated data
  vector[Nobs] y_rep;
  for (n in 1:Nobs) {
    y_rep[n] = normal_rng(mu[n], sigma);
  }
  
  // posterior expected temperature response per species
  matrix[Nspecies, N_new] f_new;
  for (m in 1:Nspecies) {
    for (n in 1:N_new) {
      f_new[m, n] = linear(x_new[n], Tmax[m], beta[m]);
    }
  }
  
  // posterior expected temperature response new species 
  vector[N_new] f_new_spec;
  real Tmax_new;
  real beta_new;
  Tmax_new = normal_rng(mu_Tmax,sigma_Tmax);
  beta_new = lognormal_rng(mu_beta,sigma_beta);
  for (n in 1:N_new) {
      f_new_spec[n] = linear(x_new[n], Tmax_new, beta_new);
  }
}
