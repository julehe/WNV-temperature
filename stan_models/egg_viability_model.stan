functions {
  // a modified quadratic function
  real quad(real x, real q, real Tmin, real Tmax) {
    real out;
    if (x>Tmin && x<Tmax) {
      out = (q/1000 * (x - Tmin) * (Tmax - x));
    } else {
      out = 0;
    }
    if (out>1) {
      out = 1;
    }
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
  real mu_Tmin;
  real mu_Tmax;
  real mu_q;
  
  // between-species standard deviations
  real<lower=0> sigma_Tmin;
  real<lower=0> sigma_Tmax;
  real<lower=0> sigma_q;
  
  // non-centered species-level parameters
  vector[Nspecies] Tmin_raw;
  vector[Nspecies] Tmax_raw;
  vector[Nspecies] q_raw;
  
  // between-experiment standard deviations
  real<lower=0> sigma_exp_Tmin;
  real<lower=0> sigma_exp_Tmax;
  real<lower=0> sigma_exp_q;
  // non-centered experiment-level parameters (deviations from species mean)
  vector[Nexp] exp_Tmin;
  vector[Nexp] exp_Tmax;
  vector[Nexp] exp_q;
  
  real<lower=0> sigma; // measurement error
}
transformed parameters{
  // centered species-level parameters
  vector[Nspecies] Tmin;
  vector[Nspecies] Tmax;
  vector[Nspecies] q;
  
  Tmin = Tmin_raw * sigma_Tmin + mu_Tmin;
  Tmax = Tmax_raw * sigma_Tmax + mu_Tmax;
  q = exp(q_raw * sigma_q + mu_q);
  
  // centered experiment-level parameters (deviations from species mean)
  vector[Nexp] exp_q_sc = exp(exp_q*sigma_exp_q);
  vector[Nexp] exp_Tmin_sc = exp_Tmin*sigma_exp_Tmin;
  vector[Nexp] exp_Tmax_sc =exp_Tmax*sigma_exp_Tmax;
  
  // likelihood mean for each observation
  vector[Nobs] mu; 
  for (n in 1:Nobs){
    mu[n] = quad(x[n], q[species_id[n]]*exp_q_sc[exp_id[n]], Tmin[species_id[n]]+exp_Tmin_sc[exp_id[n]], Tmax[species_id[n]]+exp_Tmax_sc[exp_id[n]]);
  }
}

model {
  // hyperpriors
  mu_Tmin ~ normal(10,10);
  mu_Tmax ~ normal(38,10);
  mu_q ~ normal(1.5,1);
  
  sigma_Tmin ~ gamma(7.96,2.85);
  sigma_Tmax ~ gamma(4.86,1.88);
  sigma_q ~ gamma(1.2,16.37);
  
  sigma_exp_Tmin ~ gamma(3.45,2.35);
  sigma_exp_Tmax ~ gamma(18.93,10.07);
  sigma_exp_q ~ gamma(2.79,15.63);
  
  // priors
  Tmin_raw ~ std_normal();
  Tmax_raw ~ std_normal();
  q_raw ~ std_normal();
  
  exp_Tmin ~ std_normal();
  exp_Tmax ~ std_normal();
  exp_q ~ std_normal();
  
  sigma ~ std_normal();
  
  // likelihood
  y ~ normal(mu, sigma);
}
generated quantities {
  // optimal temperature at species-level
  vector[Nspecies] Topt;
  for (n in 1:Nspecies){
    Topt[n] = 0.5 * (Tmin[n] + Tmax[n]);
  }

  // posterior predicition replicated data
  vector[Nobs] y_rep;
  for (n in 1:Nobs) {
    y_rep[n] = normal_rng(mu[n], sigma);
  }
  
  // posterior expected temperature response per species
  matrix[Nspecies, N_new] f_new;
  for (m in 1:Nspecies) {
    for (n in 1:N_new) {
      f_new[m, n] = quad(x_new[n], q[m], Tmin[m], Tmax[m]);
    }
  }
  
  // posterior expected temperature response new species 
  vector[N_new] f_new_spec;
  real Tmax_new;
  real Tmin_new;
  real q_new;
  Tmax_new = normal_rng(mu_Tmax,sigma_Tmax);
  Tmin_new = normal_rng(mu_Tmin,sigma_Tmin);
  q_new = lognormal_rng(mu_q,sigma_q);
  for (n in 1:N_new) {
      f_new_spec[n] = quad(x_new[n], q_new, Tmin_new, Tmax_new);
  }
}
