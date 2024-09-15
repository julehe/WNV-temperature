functions {
  // a linear function
  real linear(real x, real alpha, real beta) {
    real out;
    out = beta*x + alpha;
    return out;
  }
}
data {
  int<lower=0> Nobs; // number of observations
  int<lower=0> Nexp; // number of experiments
  int<lower=0> ns[Nobs]; // sample size of each observation
  vector[Nobs] x; // predictor (temperature)
  int<lower=0> y[Nobs]; // response (trait)
  
  array[Nobs] int<lower=1, upper=Nexp> exp_id; // encodes to which experiment each observation belongs
  
  int<lower=0> N_new; // number of temp. points for generated quantities
  vector[N_new] x_new; // temp. points for generated quantities
}
parameters {
  // population-level means
  real mu_alpha;
  real mu_beta;
  
  // between-experiment standard deviations
  vector[Nexp] exp_alpha;
  vector[Nexp] exp_beta;
}
transformed parameters{
  // experiment-level parameters
  vector[Nexp] alpha;
  vector[Nexp] beta;
  for (n in 1:Nexp) {
    alpha[n] = mu_alpha + exp_alpha[n];
    beta[n] = mu_beta + exp_beta[n];
  }
  
  // success probability for each observation
  vector[Nobs] mu;
  for (n in 1:Nobs) {
    mu[n] = inv_logit(linear(x[n], mu_alpha + exp_alpha[exp_id[n]], mu_beta + exp_beta[exp_id[n]]));
  }
}
model {
  // hyperpriors
  mu_beta ~ normal(0,1);
  mu_alpha ~ normal(0,10);
  
  // priors
  exp_alpha ~ normal(0,0.45);
  exp_beta ~ normal(0,0.03);
  
  // likelihood
  y ~ binomial(ns, mu);
}
generated quantities {
  // posterior predicition replicated data
  vector[Nobs] y_rep;
  for (n in 1:Nobs) {
    y_rep[n] = binomial_rng(ns[n], mu[n]);
  }
  
  // posterior expected temperature response per experiment
  vector[N_new] f_new_Dohm;
  vector[N_new] f_new_KPWN02;
  vector[N_new] f_new_KPNY99;
  for (n in 1:N_new) {
    f_new_Dohm[n] = inv_logit(linear(x_new[n], mu_alpha + exp_alpha[1], mu_beta + exp_beta[1]));
    f_new_KPWN02[n] = inv_logit(linear(x_new[n], mu_alpha + exp_alpha[2], mu_beta + exp_beta[2]));
    f_new_KPNY99[n] = inv_logit(linear(x_new[n], mu_alpha + exp_alpha[3], mu_beta + exp_beta[3]));
  }
  
  // posterior expected temperature response new experiment
  vector[N_new] f_new_spec;
  real alpha_new;
  real beta_new;
  alpha_new = normal_rng(mu_alpha,0.3);
  beta_new = normal_rng(mu_beta,0.02);
  for (n in 1:N_new) {
    f_new_spec[n] = inv_logit(linear(x_new[n], alpha_new, beta_new));
  }
}
