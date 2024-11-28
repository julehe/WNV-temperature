### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves uncertainty and sensitivity analysis.
###
### The script includes the necessary code to reproduce Figure SI11.1-11
###

# load libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(ggpubr)
library(sensitivity)

## load fitted models, extract trait temperature response samples for each species
## and population-level estimates (which we use to substitute the temperature 
## response of a species when it lacks data for a specific trait), and remove models
## afterwards to save local memory

# EIP
EIP_mean <- 28
EIP_fit <- readRDS("model_fits/EIP_fit.rds")
EIP_fit_pop <- rstan::extract(EIP_fit, permuted=T)$f_new_spec
EIP_fit_Cpip <- rstan::extract(EIP_fit, permuted=T)$f_new[,1,]
EIP_fit_Ctar <- rstan::extract(EIP_fit, permuted=T)$f_new[,4,]
rm(EIP_fit)
EIP_fit_increased <- readRDS("model_fits/EIP_fit_increased_experiment_variability.rds")
EIP_fit_pop_inc <- rstan::extract(EIP_fit_increased, permuted=T)$f_new_spec
EIP_fit_Cpip_inc <- rstan::extract(EIP_fit_increased, permuted=T)$f_new[,1,]
EIP_fit_Ctar_inc <- rstan::extract(EIP_fit_increased, permuted=T)$f_new[,4,]
rm(EIP_fit_increased)
EIP_fit_reduced <- readRDS("model_fits/EIP_fit_reduced_experiment_variability.rds")
EIP_fit_pop_red <- rstan::extract(EIP_fit_reduced, permuted=T)$f_new_spec
EIP_fit_Cpip_red <- rstan::extract(EIP_fit_reduced, permuted=T)$f_new[,1,]
EIP_fit_Ctar_red <- rstan::extract(EIP_fit_reduced, permuted=T)$f_new[,4,]
rm(EIP_fit_reduced)

# biting rate
biting_data <- read.csv("data/biting_rate_data.csv")
biting_data_mean <- mean(1/biting_data$trait)
biting_fit <- readRDS("model_fits/biting_rate_fit.rds")
biting_fit_Cpal <- rstan::extract(biting_fit, permuted=T)$f_new[,2,]
biting_fit_Cpip <- rstan::extract(biting_fit, permuted=T)$f_new[,3,]
biting_fit_Cqui <- rstan::extract(biting_fit, permuted=T)$f_new[,4,]
biting_fit_Ctar <- rstan::extract(biting_fit, permuted=T)$f_new[,5,]
biting_fit_pop <- rstan::extract(biting_fit, permuted=T)$f_new_spec
rm(biting_fit)

# Adult mosquito lifespan
lf_data <- read.csv("data/lifespan_data.csv")
lf_data_mean <- mean(lf_data$trait)
lf_fit <- readRDS("model_fits/lifespan_fit.rds")
lf_fit_Cmol <- rstan::extract(lf_fit, permuted=T)$f_new[,2,]
lf_fit_Cpal <- rstan::extract(lf_fit, permuted=T)$f_new[,3,]
lf_fit_Cpip <- rstan::extract(lf_fit, permuted=T)$f_new[,4,]
lf_fit_Cqui <- rstan::extract(lf_fit, permuted=T)$f_new[,5,]
lf_fit_Cres <- rstan::extract(lf_fit, permuted=T)$f_new[,6,]
lf_fit_Ctar <- rstan::extract(lf_fit, permuted=T)$f_new[,7,]
rm(lf_fit)

# introduce cut off at lowest observed temperature (14°C)
steps = 0.1
for(i in 1:(14/steps)){
  lf_fit_Cmol[,i] <- lf_fit_Cmol[,(14/steps+1)]
  lf_fit_Cpal[,i] <- lf_fit_Cpal[,(14/steps+1)]
  lf_fit_Cpip[,i] <- lf_fit_Cpip[,(14/steps+1)]
  lf_fit_Cqui[,i] <- lf_fit_Cqui[,(14/steps+1)]
  lf_fit_Cres[,i] <- lf_fit_Cres[,(14/steps+1)]
  lf_fit_Ctar[,i] <- lf_fit_Ctar[,(14/steps+1)]
}

# Egg development rate
dev_egg_data <- read.csv("data/egg_development_data.csv")
dev_egg_data_mean <- mean(1/dev_egg_data$trait)
dev_egg_fit <- readRDS("model_fits/egg_development_fit.rds")
dev_egg_fit_Cmol <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,2,]
dev_egg_fit_Cpal <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,3,]
dev_egg_fit_Cpip <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,4,]
dev_egg_fit_Cqui <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,5,]
dev_egg_fit_Cres <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,6,]
rm(dev_egg_fit)

# Juvenile development rate
dev_data <- read.csv("data/juvenile_development_data.csv")
dev_data_mean <- mean(1/dev_data$trait)
dev_fit <- readRDS("model_fits/juvenile_development_fit.rds")
dev_fit_Cmol <- rstan::extract(dev_fit, permuted=T)$f_new[,5,]
dev_fit_Cpal <- rstan::extract(dev_fit, permuted=T)$f_new[,6,]
dev_fit_Cpip <- rstan::extract(dev_fit, permuted=T)$f_new[,7,]
dev_fit_Cqui <- rstan::extract(dev_fit, permuted=T)$f_new[,8,]
dev_fit_Cres <- rstan::extract(dev_fit, permuted=T)$f_new[,9,]
dev_fit_Ctar <- rstan::extract(dev_fit, permuted=T)$f_new[,11,]
rm(dev_fit)

# Egg viability
egg_viability_data <- read.csv("data/egg_viability_data.csv")
egg_viability_data_mean <- mean(egg_viability_data$trait)
egg_viability_fit <- readRDS("model_fits/egg_viability_fit.rds")
egg_viability_fit_Cmol <- rstan::extract(egg_viability_fit, permuted=T)$f_new[,1,]
egg_viability_fit_Cpal <- rstan::extract(egg_viability_fit, permuted=T)$f_new[,2,]
egg_viability_fit_Cqui <- rstan::extract(egg_viability_fit, permuted=T)$f_new[,3,]
egg_viability_fit_pop <- rstan::extract(egg_viability_fit, permuted=T)$f_new_spec
rm(egg_viability_fit)

# Mosquito infection probability
infprob_data <- read.csv("data/infection_probability_data.csv")
infprob_data_mean <- mean(infprob_data$trait)
infprob_fit <- readRDS("model_fits/infection_probability_fit.rds")
infprob_fit_pop <- rstan::extract(infprob_fit, permuted=T)$f_new_spec
rm(infprob_fit)

# Juvenile survival
sur_data <- read.csv("data/juvenile_survival_data.csv")
sur_data_mean <- mean(sur_data$trait)
sur_fit <- readRDS("model_fits/juvenile_survival_fit.rds")
sur_fit_Cmol <- rstan::extract(sur_fit, permuted=T)$f_new[,7,]
sur_fit_Cpal <- rstan::extract(sur_fit, permuted=T)$f_new[,8,]
sur_fit_Cpip <- rstan::extract(sur_fit, permuted=T)$f_new[,9,]
sur_fit_Cqui <- rstan::extract(sur_fit, permuted=T)$f_new[,10,]
sur_fit_Cres <- rstan::extract(sur_fit, permuted=T)$f_new[,11,]
sur_fit_Ctar <- rstan::extract(sur_fit, permuted=T)$f_new[,13,]
rm(sur_fit)

set.seed(123)
X_index1 <- sample(x = 1:8000, size = 4000, replace = FALSE)
X_index2 <- setdiff(1:8000, X_index1)
# set temperature points at which model outputs were generated on
temp = seq(0,45,steps)
# small constant used to avoid division by zero
cc <- 0.0000001
# calculate mean number of eggs per raft
ER_data <- read.csv("data/eggs_per_raft_data.csv")
ER = round(mean(ER_data$trait))
# assumed sex ratio at adult emergence
omega = 0.5

## some helper functions

# function evaluating the mosquito abundance approximation used in the main text
M_main_f <- function(omega, surJ, ER, a, EV, lf, devJ){
  EFD = ER * a
  surEJ = surJ * EV
  ifelse(1/(lf * omega * EFD * surEJ + cc) <1,
         (omega^2 * EFD * EV * devJ^2 * lf^2) * (1 - 1/(lf * omega * EFD * surEJ + cc)),
         0)
}

# function evaluating relative R0 using the main mosquito abundance approximation
R0_main_f = function(a, b, EIP, lf, omega, surJ, ER, EV, devJ){
  muM = 1/(lf+cc)
  (a^2 * b * exp(-EIP*muM) * M_main_f(omega, surJ, ER, a, EV, lf, devJ))/(muM)
}
R0_main_sobol = function(X){
  a = X[,1]
  b = X[,2]
  EIP = X[,3]
  lf = X[,4]
  surJ = X[,5]
  EV = X[,6]
  devJ = X[,7]
  muM = 1/(lf+cc)
  (a^2 * b * exp(-EIP*muM) * M_main_f(omega, surJ, ER, a, EV, lf, devJ))/(muM)
}

# function to calculate Tmin (lower temperature where R0 becomes zero)
Tmin  = function(x){
  index_list = which(x>0) 
  temp[index_list[1] - 1]
}

# function to calculate Tmax (upper temperature where R0 becomes zero)
Tmax  = function(x){
  index_list = which(x>0)
  temp[index_list[length(index_list)] + 1]
}

## relative R0 calculations for each species 

# Cx. pipiens

# calculate R0 samples from the trait samples
R0_main_Cpip = R0_main_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, 
                         lf_fit_Cpip, omega, sur_fit_Cpip, ER, 
                         egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip <= 0) == ncol(R0_main_Cpip))){
  print("only zeros")
  R0_main_Cpip <- R0_main_Cpip[-which(rowSums(R0_main_Cpip <= 0) == ncol(R0_main_Cpip)),]  
}

# calculate statistics of relative R0 at each temperature
R0_main_Cpip_mean = apply(R0_main_Cpip, 2, mean)
R0_main_Cpip_0025 = apply(R0_main_Cpip, 2, quantile, probs=c(0.025))
R0_main_Cpip_0975 = apply(R0_main_Cpip, 2, quantile, probs=c(0.975))

# calculate Topt for each sample and calculate statistics
R0_main_Cpip_peaks = sapply(apply(R0_main_Cpip,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_peaks_mean = mean(R0_main_Cpip_peaks)
R0_main_Cpip_peaks_median = median(R0_main_Cpip_peaks)
R0_main_Cpip_peaks_0025 = quantile(R0_main_Cpip_peaks, probs=c(0.025))
R0_main_Cpip_peaks_0975 = quantile(R0_main_Cpip_peaks, probs=c(0.975))

# calculate Tmin for each sample and calculate statistics
R0_main_Cpip_Tmin = apply(R0_main_Cpip,1,FUN=Tmin)
R0_main_Cpip_Tmin_mean = mean(R0_main_Cpip_Tmin)
R0_main_Cpip_Tmin_median = median(R0_main_Cpip_Tmin)
R0_main_Cpip_Tmin_0025 = quantile(R0_main_Cpip_Tmin, probs=c(0.025))
R0_main_Cpip_Tmin_0975 = quantile(R0_main_Cpip_Tmin, probs=c(0.975))

# calculate Tmax for each sample and calculate statistics
R0_main_Cpip_Tmax = apply(R0_main_Cpip,1,FUN=Tmax)
R0_main_Cpip_Tmax_mean = mean(R0_main_Cpip_Tmax)
R0_main_Cpip_Tmax_median = median(R0_main_Cpip_Tmax)
R0_main_Cpip_Tmax_0025 = quantile(R0_main_Cpip_Tmax, probs=c(0.025))
R0_main_Cpip_Tmax_0975 = quantile(R0_main_Cpip_Tmax, probs=c(0.975))

# calculate R0 samples again but leaving one trait constant

# Biting rate

R0_main_Cpip_biting_const = R0_main_f(biting_data_mean, infprob_fit_pop, EIP_fit_Cpip, 
                                      lf_fit_Cpip, omega, sur_fit_Cpip, ER, 
                                      egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_biting_const <= 0) == ncol(R0_main_Cpip_biting_const))){
  print("only zeros")
  R0_main_Cpip_biting_const <- R0_main_Cpip[-which(rowSums(R0_main_Cpip_biting_const <= 0) == ncol(R0_main_Cpip_biting_const)),]  
}

R0_main_Cpip_biting_const_mean = apply(R0_main_Cpip_biting_const, 2, mean)

R0_main_Cpip_biting_const_peaks = sapply(apply(R0_main_Cpip_biting_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_biting_const_peaks_mean = mean(R0_main_Cpip_biting_const_peaks)
R0_main_Cpip_biting_const_peaks_0025 = quantile(R0_main_Cpip_biting_const_peaks, probs=c(0.025))
R0_main_Cpip_biting_const_peaks_0975 = quantile(R0_main_Cpip_biting_const_peaks, probs=c(0.975))

R0_main_Cpip_biting_const_Tmax = apply(R0_main_Cpip_biting_const,1,FUN=Tmax)
R0_main_Cpip_biting_const_Tmax_mean = mean(R0_main_Cpip_biting_const_Tmax)
R0_main_Cpip_biting_const_Tmax_0025 = quantile(R0_main_Cpip_biting_const_Tmax, probs=c(0.025))
R0_main_Cpip_biting_const_Tmax_0975 = quantile(R0_main_Cpip_biting_const_Tmax, probs=c(0.975))

R0_main_Cpip_biting_const_Tmin = apply(R0_main_Cpip_biting_const,1,FUN=Tmin)
R0_main_Cpip_biting_const_Tmin_mean = mean(R0_main_Cpip_biting_const_Tmin)
R0_main_Cpip_biting_const_Tmin_0025 = quantile(R0_main_Cpip_biting_const_Tmin, probs=c(0.025))
R0_main_Cpip_biting_const_Tmin_0975 = quantile(R0_main_Cpip_biting_const_Tmin, probs=c(0.975))

# Infection probability

R0_main_Cpip_infprob_const = R0_main_f(biting_fit_Cpip, infprob_data_mean, EIP_fit_Cpip, 
                                       lf_fit_Cpip, omega, sur_fit_Cpip, ER, 
                                       egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_infprob_const <= 0) == ncol(R0_main_Cpip_infprob_const))){
  print("only zeros")
  R0_main_Cpip_infprob_const <- R0_main_Cpip[-which(rowSums(R0_main_Cpip_infprob_const <= 0) == ncol(R0_main_Cpip_infprob_const)),]  
}

R0_main_Cpip_infprob_const_mean = apply(R0_main_Cpip_infprob_const, 2, mean)

R0_main_Cpip_infprob_const_peaks = sapply(apply(R0_main_Cpip_infprob_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_infprob_const_peaks_mean = mean(R0_main_Cpip_infprob_const_peaks)
R0_main_Cpip_infprob_const_peaks_0025 = quantile(R0_main_Cpip_infprob_const_peaks, probs=c(0.025))
R0_main_Cpip_infprob_const_peaks_0975 = quantile(R0_main_Cpip_infprob_const_peaks, probs=c(0.975))

R0_main_Cpip_infprob_const_Tmax = apply(R0_main_Cpip_infprob_const,1,FUN=Tmax)
R0_main_Cpip_infprob_const_Tmax_mean = mean(R0_main_Cpip_infprob_const_Tmax)
R0_main_Cpip_infprob_const_Tmax_0025 = quantile(R0_main_Cpip_infprob_const_Tmax, probs=c(0.025))
R0_main_Cpip_infprob_const_Tmax_0975 = quantile(R0_main_Cpip_infprob_const_Tmax, probs=c(0.975))

R0_main_Cpip_infprob_const_Tmin = apply(R0_main_Cpip_infprob_const,1,FUN=Tmin)
R0_main_Cpip_infprob_const_Tmin_mean = mean(R0_main_Cpip_infprob_const_Tmin)
R0_main_Cpip_infprob_const_Tmin_0025 = quantile(R0_main_Cpip_infprob_const_Tmin, probs=c(0.025))
R0_main_Cpip_infprob_const_Tmin_0975 = quantile(R0_main_Cpip_infprob_const_Tmin, probs=c(0.975))

# EIP

R0_main_Cpip_EIP_const = R0_main_f(biting_fit_Cpip, infprob_fit_pop, EIP_mean, 
                                   lf_fit_Cpip, omega, sur_fit_Cpip, ER, 
                                   egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_EIP_const <= 0) == ncol(R0_main_Cpip_EIP_const))){
  print("only zeros")
  R0_main_Cpip_EIP_const <- R0_main_Cpip[-which(rowSums(R0_main_Cpip_EIP_const <= 0) == ncol(R0_main_Cpip_EIP_const)),]  
}

R0_main_Cpip_EIP_const_mean = apply(R0_main_Cpip_EIP_const, 2, mean)

R0_main_Cpip_EIP_const_peaks = sapply(apply(R0_main_Cpip_EIP_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_EIP_const_peaks_mean = mean(R0_main_Cpip_EIP_const_peaks)
R0_main_Cpip_EIP_const_peaks_0025 = quantile(R0_main_Cpip_EIP_const_peaks, probs=c(0.025))
R0_main_Cpip_EIP_const_peaks_0975 = quantile(R0_main_Cpip_EIP_const_peaks, probs=c(0.975))

R0_main_Cpip_EIP_const_Tmax = apply(R0_main_Cpip_EIP_const,1,FUN=Tmax)
R0_main_Cpip_EIP_const_Tmax_mean = mean(R0_main_Cpip_EIP_const_Tmax)
R0_main_Cpip_EIP_const_Tmax_0025 = quantile(R0_main_Cpip_EIP_const_Tmax, probs=c(0.025))
R0_main_Cpip_EIP_const_Tmax_0975 = quantile(R0_main_Cpip_EIP_const_Tmax, probs=c(0.975))

R0_main_Cpip_EIP_const_Tmin = apply(R0_main_Cpip_EIP_const,1,FUN=Tmin)
R0_main_Cpip_EIP_const_Tmin_mean = mean(R0_main_Cpip_EIP_const_Tmin)
R0_main_Cpip_EIP_const_Tmin_0025 = quantile(R0_main_Cpip_EIP_const_Tmin, probs=c(0.025))
R0_main_Cpip_EIP_const_Tmin_0975 = quantile(R0_main_Cpip_EIP_const_Tmin, probs=c(0.975))

# Adult lifespan

R0_main_Cpip_lf_const = R0_main_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, 
                                  lf_data_mean, omega, sur_fit_Cpip, ER, 
                                  egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_lf_const <= 0) == ncol(R0_main_Cpip_lf_const))){
  print("only zeros")
  R0_main_Cpip_lf_const <- R0_main_Cpip[-which(rowSums(R0_main_Cpip_lf_const <= 0) == ncol(R0_main_Cpip_lf_const)),]  
}

R0_main_Cpip_lf_const_mean = apply(R0_main_Cpip_lf_const, 2, mean)

R0_main_Cpip_lf_const_peaks = sapply(apply(R0_main_Cpip_lf_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_lf_const_peaks_mean = mean(R0_main_Cpip_lf_const_peaks)
R0_main_Cpip_lf_const_peaks_0025 = quantile(R0_main_Cpip_lf_const_peaks, probs=c(0.025))
R0_main_Cpip_lf_const_peaks_0975 = quantile(R0_main_Cpip_lf_const_peaks, probs=c(0.975))

R0_main_Cpip_lf_const_Tmax = apply(R0_main_Cpip_lf_const,1,FUN=Tmax)
R0_main_Cpip_lf_const_Tmax_mean = mean(R0_main_Cpip_lf_const_Tmax)
R0_main_Cpip_lf_const_Tmax_0025 = quantile(R0_main_Cpip_lf_const_Tmax, probs=c(0.025))
R0_main_Cpip_lf_const_Tmax_0975 = quantile(R0_main_Cpip_lf_const_Tmax, probs=c(0.975))

R0_main_Cpip_lf_const_Tmin = apply(R0_main_Cpip_lf_const,1,FUN=Tmin)
R0_main_Cpip_lf_const_Tmin_mean = mean(R0_main_Cpip_lf_const_Tmin)
R0_main_Cpip_lf_const_Tmin_0025 = quantile(R0_main_Cpip_lf_const_Tmin, probs=c(0.025))
R0_main_Cpip_lf_const_Tmin_0975 = quantile(R0_main_Cpip_lf_const_Tmin, probs=c(0.975))

# Juvenile survival probability

R0_main_Cpip_sur_const = R0_main_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, 
                                   lf_fit_Cpip, omega, sur_data_mean, ER, 
                                   egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_sur_const <= 0) == ncol(R0_main_Cpip_sur_const))){
  print("only zeros")
  R0_main_Cpip_sur_const <- R0_main_Cpip[-which(rowSums(R0_main_Cpip_sur_const <= 0) == ncol(R0_main_Cpip_sur_const)),]  
}

R0_main_Cpip_sur_const_mean = apply(R0_main_Cpip_sur_const, 2, mean)

R0_main_Cpip_sur_const_peaks = sapply(apply(R0_main_Cpip_sur_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_sur_const_peaks_mean = mean(R0_main_Cpip_sur_const_peaks)
R0_main_Cpip_sur_const_peaks_0025 = quantile(R0_main_Cpip_sur_const_peaks, probs=c(0.025))
R0_main_Cpip_sur_const_peaks_0975 = quantile(R0_main_Cpip_sur_const_peaks, probs=c(0.975))

R0_main_Cpip_sur_const_Tmax = apply(R0_main_Cpip_sur_const,1,FUN=Tmax)
R0_main_Cpip_sur_const_Tmax_mean = mean(R0_main_Cpip_sur_const_Tmax)
R0_main_Cpip_sur_const_Tmax_0025 = quantile(R0_main_Cpip_sur_const_Tmax, probs=c(0.025))
R0_main_Cpip_sur_const_Tmax_0975 = quantile(R0_main_Cpip_sur_const_Tmax, probs=c(0.975))

R0_main_Cpip_sur_const_Tmin = apply(R0_main_Cpip_sur_const,1,FUN=Tmin)
R0_main_Cpip_sur_const_Tmin_mean = mean(R0_main_Cpip_sur_const_Tmin)
R0_main_Cpip_sur_const_Tmin_0025 = quantile(R0_main_Cpip_sur_const_Tmin, probs=c(0.025))
R0_main_Cpip_sur_const_Tmin_0975 = quantile(R0_main_Cpip_sur_const_Tmin, probs=c(0.975))

# Egg viability

R0_main_Cpip_EV_const = R0_main_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, 
                                  lf_fit_Cpip, omega, sur_fit_Cpip, ER, 
                                  egg_viability_data_mean, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_EV_const <= 0) == ncol(R0_main_Cpip_EV_const))){
  print("only zeros")
  R0_main_Cpip_EV_const <- R0_main_Cpip[-which(rowSums(R0_main_Cpip_EV_const <= 0) == ncol(R0_main_Cpip_EV_const)),]  
}

R0_main_Cpip_EV_const_mean = apply(R0_main_Cpip_EV_const, 2, mean)

R0_main_Cpip_EV_const_peaks = sapply(apply(R0_main_Cpip_EV_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_EV_const_peaks_mean = mean(R0_main_Cpip_EV_const_peaks)
R0_main_Cpip_EV_const_peaks_0025 = quantile(R0_main_Cpip_EV_const_peaks, probs=c(0.025))
R0_main_Cpip_EV_const_peaks_0975 = quantile(R0_main_Cpip_EV_const_peaks, probs=c(0.975))

R0_main_Cpip_EV_const_Tmax = apply(R0_main_Cpip_EV_const,1,FUN=Tmax)
R0_main_Cpip_EV_const_Tmax_mean = mean(R0_main_Cpip_EV_const_Tmax)
R0_main_Cpip_EV_const_Tmax_0025 = quantile(R0_main_Cpip_EV_const_Tmax, probs=c(0.025))
R0_main_Cpip_EV_const_Tmax_0975 = quantile(R0_main_Cpip_EV_const_Tmax, probs=c(0.975))

R0_main_Cpip_EV_const_Tmin = apply(R0_main_Cpip_EV_const,1,FUN=Tmin)
R0_main_Cpip_EV_const_Tmin_mean = mean(R0_main_Cpip_EV_const_Tmin)
R0_main_Cpip_EV_const_Tmin_0025 = quantile(R0_main_Cpip_EV_const_Tmin, probs=c(0.025))
R0_main_Cpip_EV_const_Tmin_0975 = quantile(R0_main_Cpip_EV_const_Tmin, probs=c(0.975))

# Juvenile development rate

R0_main_Cpip_dev_const = R0_main_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, 
                                   lf_fit_Cpip, omega, sur_fit_Cpip, ER, 
                                   egg_viability_fit_pop, dev_data_mean)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_dev_const <= 0) == ncol(R0_main_Cpip_dev_const))){
  print("only zeros")
  R0_main_Cpip_dev_const <- R0_main_Cpip[-which(rowSums(R0_main_Cpip_dev_const <= 0) == ncol(R0_main_Cpip_dev_const)),]  
}

R0_main_Cpip_dev_const_mean = apply(R0_main_Cpip_dev_const, 2, mean)

R0_main_Cpip_dev_const_peaks = sapply(apply(R0_main_Cpip_dev_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_dev_const_peaks_mean = mean(R0_main_Cpip_dev_const_peaks)
R0_main_Cpip_dev_const_peaks_0025 = quantile(R0_main_Cpip_dev_const_peaks, probs=c(0.025))
R0_main_Cpip_dev_const_peaks_0975 = quantile(R0_main_Cpip_dev_const_peaks, probs=c(0.975))

R0_main_Cpip_dev_const_Tmax = apply(R0_main_Cpip_dev_const,1,FUN=Tmax)
R0_main_Cpip_dev_const_Tmax_mean = mean(R0_main_Cpip_dev_const_Tmax)
R0_main_Cpip_dev_const_Tmax_0025 = quantile(R0_main_Cpip_dev_const_Tmax, probs=c(0.025))
R0_main_Cpip_dev_const_Tmax_0975 = quantile(R0_main_Cpip_dev_const_Tmax, probs=c(0.975))

R0_main_Cpip_dev_const_Tmin = apply(R0_main_Cpip_dev_const,1,FUN=Tmin)
R0_main_Cpip_dev_const_Tmin_mean = mean(R0_main_Cpip_dev_const_Tmin)
R0_main_Cpip_dev_const_Tmin_0025 = quantile(R0_main_Cpip_dev_const_Tmin, probs=c(0.025))
R0_main_Cpip_dev_const_Tmin_0975 = quantile(R0_main_Cpip_dev_const_Tmin, probs=c(0.975))

# Summarize results for Cx. pipiens in dataframe

df_R0_main_Cpip_const <- data.frame(x = temp,
                                    y = c(R0_main_Cpip_biting_const_mean/max(R0_main_Cpip_biting_const_mean),
                                          R0_main_Cpip_infprob_const_mean/max(R0_main_Cpip_infprob_const_mean),
                                          R0_main_Cpip_EIP_const_mean/max(R0_main_Cpip_EIP_const_mean),
                                          R0_main_Cpip_lf_const_mean/max(R0_main_Cpip_lf_const_mean),
                                          R0_main_Cpip_sur_const_mean/max(R0_main_Cpip_sur_const_mean),
                                          R0_main_Cpip_EV_const_mean/max(R0_main_Cpip_EV_const_mean),
                                          R0_main_Cpip_dev_const_mean/max(R0_main_Cpip_dev_const_mean),
                                          R0_main_Cpip_mean/max(R0_main_Cpip_mean)
                                    ),
                                    trait = c(rep("3Biting rate", length(temp)),
                                              rep("2Inf. prob.", length(temp)),
                                              rep("1EIP", length(temp)),
                                              rep("5Lifespan", length(temp)),
                                              rep("6Juv. survival", length(temp)),
                                              rep("4Egg viab.", length(temp)),
                                              rep("7Juv. dev. rate", length(temp)),
                                              rep("8None", length(temp))
                                    )
)

df_R0_main_Cpip_const_stats <- data.frame(peaks = c(R0_main_Cpip_biting_const_peaks_mean,
                                                    R0_main_Cpip_infprob_const_peaks_mean,
                                                    R0_main_Cpip_EIP_const_peaks_mean,
                                                    R0_main_Cpip_lf_const_peaks_mean,
                                                    R0_main_Cpip_sur_const_peaks_mean,
                                                    R0_main_Cpip_EV_const_peaks_mean,
                                                    R0_main_Cpip_dev_const_peaks_mean,
                                                    R0_main_Cpip_peaks_mean),
                                          peaks_lower = c(R0_main_Cpip_biting_const_peaks_0025,
                                                          R0_main_Cpip_infprob_const_peaks_0025,
                                                          R0_main_Cpip_EIP_const_peaks_0025,
                                                          R0_main_Cpip_lf_const_peaks_0025,
                                                          R0_main_Cpip_sur_const_peaks_0025,
                                                          R0_main_Cpip_EV_const_peaks_0025,
                                                          R0_main_Cpip_dev_const_peaks_0025,
                                                          R0_main_Cpip_peaks_0025),
                                          peaks_upper = c(R0_main_Cpip_biting_const_peaks_0975,
                                                          R0_main_Cpip_infprob_const_peaks_0975,
                                                          R0_main_Cpip_EIP_const_peaks_0975,
                                                          R0_main_Cpip_lf_const_peaks_0975,
                                                          R0_main_Cpip_sur_const_peaks_0975,
                                                          R0_main_Cpip_EV_const_peaks_0975,
                                                          R0_main_Cpip_dev_const_peaks_0975,
                                                          R0_main_Cpip_peaks_0975),
                                          Tmin = c(R0_main_Cpip_biting_const_Tmin_mean,
                                                   R0_main_Cpip_infprob_const_Tmin_mean,
                                                   R0_main_Cpip_EIP_const_Tmin_mean,
                                                   R0_main_Cpip_lf_const_Tmin_mean,
                                                   R0_main_Cpip_sur_const_Tmin_mean,
                                                   R0_main_Cpip_EV_const_Tmin_mean,
                                                   R0_main_Cpip_dev_const_Tmin_mean,
                                                   R0_main_Cpip_Tmin_mean),
                                          Tmin_lower = c(R0_main_Cpip_biting_const_Tmin_0025,
                                                         R0_main_Cpip_infprob_const_Tmin_0025,
                                                         R0_main_Cpip_EIP_const_Tmin_0025,
                                                         R0_main_Cpip_lf_const_Tmin_0025,
                                                         R0_main_Cpip_sur_const_Tmin_0025,
                                                         R0_main_Cpip_EV_const_Tmin_0025,
                                                         R0_main_Cpip_dev_const_Tmin_0025,
                                                         R0_main_Cpip_Tmin_0025),
                                          Tmin_upper = c(R0_main_Cpip_biting_const_Tmin_0975,
                                                         R0_main_Cpip_infprob_const_Tmin_0975,
                                                         R0_main_Cpip_EIP_const_Tmin_0975,
                                                         R0_main_Cpip_lf_const_Tmin_0975,
                                                         R0_main_Cpip_sur_const_Tmin_0975,
                                                         R0_main_Cpip_EV_const_Tmin_0975,
                                                         R0_main_Cpip_dev_const_Tmin_0975,
                                                         R0_main_Cpip_Tmin_0975),
                                          Tmax = c(R0_main_Cpip_biting_const_Tmax_mean,
                                                   R0_main_Cpip_infprob_const_Tmax_mean,
                                                   R0_main_Cpip_EIP_const_Tmax_mean,
                                                   R0_main_Cpip_lf_const_Tmax_mean,
                                                   R0_main_Cpip_sur_const_Tmax_mean,
                                                   R0_main_Cpip_EV_const_Tmax_mean,
                                                   R0_main_Cpip_dev_const_Tmax_mean,
                                                   R0_main_Cpip_Tmax_mean),
                                          Tmax_lower = c(R0_main_Cpip_biting_const_Tmax_0025,
                                                         R0_main_Cpip_infprob_const_Tmax_0025,
                                                         R0_main_Cpip_EIP_const_Tmax_0025,
                                                         R0_main_Cpip_lf_const_Tmax_0025,
                                                         R0_main_Cpip_sur_const_Tmax_0025,
                                                         R0_main_Cpip_EV_const_Tmax_0025,
                                                         R0_main_Cpip_dev_const_Tmax_0025,
                                                         R0_main_Cpip_Tmax_0025),
                                          Tmax_upper = c(R0_main_Cpip_biting_const_Tmax_0975,
                                                         R0_main_Cpip_infprob_const_Tmax_0975,
                                                         R0_main_Cpip_EIP_const_Tmax_0975,
                                                         R0_main_Cpip_lf_const_Tmax_0975,
                                                         R0_main_Cpip_sur_const_Tmax_0975,
                                                         R0_main_Cpip_EV_const_Tmax_0975,
                                                         R0_main_Cpip_dev_const_Tmax_0975,
                                                         R0_main_Cpip_Tmax_0975),
                                          trait = c("3Biting rate", "2Inf. prob.", "1EIP",
                                                    "5Lifespan", "6Juv. survival", "4Egg viab.", "7Juv. dev. rate",
                                                    "8None")
)

plot_Cpip_const <- ggplot() +
  geom_line(df_R0_main_Cpip_const, mapping = aes(x = x, y = y, color = trait), linewidth=0.6) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. pipiens")))) +
  labs(x = "Temperature (°C)",
       color = "Constant trait", 
       title = expression("Mean temperature response of " * R[0]^rel)) +
  scale_color_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                                  "Lifespan", "Juv. survival", "Juv. dev. rate",
                                  "None")) +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_Cpip_stats <- ggplot(df_R0_main_Cpip_const_stats) +
  geom_pointrange(aes(x=peaks, xmin=peaks_lower, xmax=peaks_upper, y=trait, color=trait), linewidth=1.2, alpha=0.5, orientation = "y") +
  geom_pointrange(aes(x=Tmin, xmin=Tmin_lower, xmax=Tmin_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  geom_pointrange(aes(x=Tmax, xmin=Tmax_lower, xmax=Tmax_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  guides(color = "none") +
  theme_bw() +
  labs(x = "Temperature (°C)",
       title = "      Temperature limits and optimal temperature") +
  #scale_y_discrete(expand = c(0.1, 0.1)) +
  scale_y_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                              "Lifespan", "Juv. survival", "Juv. dev. rate",
                              "None")) +
  theme(plot.margin = unit(c(0.4, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cpip_const, plot_Cpip_stats, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. pipiens", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/R0_const_Cpip.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Sobol' analysis

r = 151:301
n = length(r)

sobols_first_Cpip <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

sobols_total_Cpip <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

for(i in 1:n){
  input_matrix <- data.frame(a = biting_fit_Cpip[,r[i]], bM = infprob_fit_pop[,r[i]], EIP = EIP_fit_Cpip[,r[i]], 
                             lf = lf_fit_Cpip[,r[i]], surJ = sur_fit_Cpip[,r[i]], 
                             EV= egg_viability_fit_pop[,r[i]], devJ = dev_fit_Cpip[,r[i]])
  
  X1 <- input_matrix[X_index1, ]
  X2 <- input_matrix[X_index2, ]
  
  sobol_result <- sobolmartinez(model = R0_main_sobol, 
                                X1 = X1, 
                                X2 = X2, 
                                nboot = 100)
  
  sobols_first_Cpip[i,-1] <- sobol_result$S$original
  sobols_total_Cpip[i,-1] <- sobol_result$T$original
}
sobols_first_Cpip[sobols_first_Cpip < 0] <- 0
sobols_total_Cpip[sobols_total_Cpip < 0] <- 0

sobols_first_long_Cpip <- sobols_first_Cpip %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cpip_sobol_first <- ggplot(sobols_first_long_Cpip) +
  geom_line(mapping = aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) + 
  labs(x = "Temperature (°C)",
       title = "First-order Sobol' indices",
       color = "Trait") +
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

sobols_total_long_Cpip <- sobols_total_Cpip %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cpip_sobol_total <-ggplot(sobols_total_long_Cpip) +
  geom_line(aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) +  
  labs(x = "Temperature (°C)",
       title = "Total-effect Sobol' indices",
       color = "Trait") + 
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cpip_sobol_first, plot_Cpip_sobol_total, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. pipiens", 
                                        size = 11, face = "bold.italic"))
plot

#ggsave("Figures/Sobol_Cpip.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")


# Cx. quin.

# calculate R0 samples from the trait samples
R0_main_Cqui = R0_main_f(biting_fit_Cqui, infprob_fit_pop, EIP_fit_pop, 
                         lf_fit_Cqui, omega, sur_fit_Cqui, ER, 
                         egg_viability_fit_Cqui, dev_fit_Cqui)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cqui <= 0) == ncol(R0_main_Cqui))){
  print("only zeros")
  R0_main_Cqui <- R0_main_Cqui[-which(rowSums(R0_main_Cqui <= 0) == ncol(R0_main_Cqui)),]  
}

# calculate statistics of relative R0 at each temperature
R0_main_Cqui_mean = apply(R0_main_Cqui, 2, mean)
R0_main_Cqui_0025 = apply(R0_main_Cqui, 2, quantile, probs=c(0.025))
R0_main_Cqui_0975 = apply(R0_main_Cqui, 2, quantile, probs=c(0.975))

# calculate Topt for each sample and calculate statistics
R0_main_Cqui_peaks = sapply(apply(R0_main_Cqui,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_peaks_mean = mean(R0_main_Cqui_peaks)
R0_main_Cqui_peaks_median = median(R0_main_Cqui_peaks)
R0_main_Cqui_peaks_0025 = quantile(R0_main_Cqui_peaks, probs=c(0.025))
R0_main_Cqui_peaks_0975 = quantile(R0_main_Cqui_peaks, probs=c(0.975))

# calculate Tmin for each sample and calculate statistics
R0_main_Cqui_Tmin = apply(R0_main_Cqui,1,FUN=Tmin)
R0_main_Cqui_Tmin_mean = mean(R0_main_Cqui_Tmin)
R0_main_Cqui_Tmin_median = median(R0_main_Cqui_Tmin)
R0_main_Cqui_Tmin_0025 = quantile(R0_main_Cqui_Tmin, probs=c(0.025))
R0_main_Cqui_Tmin_0975 = quantile(R0_main_Cqui_Tmin, probs=c(0.975))

# calculate Tmax for each sample and calculate statistics
R0_main_Cqui_Tmax = apply(R0_main_Cqui,1,FUN=Tmax)
R0_main_Cqui_Tmax_mean = mean(R0_main_Cqui_Tmax)
R0_main_Cqui_Tmax_median = median(R0_main_Cqui_Tmax)
R0_main_Cqui_Tmax_0025 = quantile(R0_main_Cqui_Tmax, probs=c(0.025))
R0_main_Cqui_Tmax_0975 = quantile(R0_main_Cqui_Tmax, probs=c(0.975))

# calculate R0 samples again but leaving one trait constant

# Biting rate

R0_main_Cqui_biting_const = R0_main_f(biting_data_mean, infprob_fit_pop, EIP_fit_pop, 
                                      lf_fit_Cqui, omega, sur_fit_Cqui, ER, 
                                      egg_viability_fit_Cqui, dev_fit_Cqui)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cqui_biting_const <= 0) == ncol(R0_main_Cqui_biting_const))){
  print("only zeros")
  R0_main_Cqui_biting_const <- R0_main_Cqui_biting_const[-which(rowSums(R0_main_Cqui_biting_const <= 0) == ncol(R0_main_Cqui_biting_const)),]  
}

R0_main_Cqui_biting_const_mean = apply(R0_main_Cqui_biting_const, 2, mean)

R0_main_Cqui_biting_const_peaks = sapply(apply(R0_main_Cqui_biting_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_biting_const_peaks_mean = mean(R0_main_Cqui_biting_const_peaks)
R0_main_Cqui_biting_const_peaks_0025 = quantile(R0_main_Cqui_biting_const_peaks, probs=c(0.025))
R0_main_Cqui_biting_const_peaks_0975 = quantile(R0_main_Cqui_biting_const_peaks, probs=c(0.975))

R0_main_Cqui_biting_const_Tmax = apply(R0_main_Cqui_biting_const,1,FUN=Tmax)
R0_main_Cqui_biting_const_Tmax_mean = mean(R0_main_Cqui_biting_const_Tmax)
R0_main_Cqui_biting_const_Tmax_0025 = quantile(R0_main_Cqui_biting_const_Tmax, probs=c(0.025))
R0_main_Cqui_biting_const_Tmax_0975 = quantile(R0_main_Cqui_biting_const_Tmax, probs=c(0.975))

R0_main_Cqui_biting_const_Tmin = apply(R0_main_Cqui_biting_const,1,FUN=Tmin)
R0_main_Cqui_biting_const_Tmin_mean = mean(R0_main_Cqui_biting_const_Tmin)
R0_main_Cqui_biting_const_Tmin_0025 = quantile(R0_main_Cqui_biting_const_Tmin, probs=c(0.025))
R0_main_Cqui_biting_const_Tmin_0975 = quantile(R0_main_Cqui_biting_const_Tmin, probs=c(0.975))

# Infection probability

R0_main_Cqui_infprob_const = R0_main_f(biting_fit_Cqui, infprob_data_mean, EIP_fit_pop, 
                                       lf_fit_Cqui, omega, sur_fit_Cqui, ER, 
                                       egg_viability_fit_Cqui, dev_fit_Cqui)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cqui_infprob_const <= 0) == ncol(R0_main_Cqui_infprob_const))){
  print("only zeros")
  R0_main_Cqui_infprob_const <- R0_main_Cqui_infprob_const[-which(rowSums(R0_main_Cqui_infprob_const <= 0) == ncol(R0_main_Cqui_infprob_const)),]  
}

R0_main_Cqui_infprob_const_mean = apply(R0_main_Cqui_infprob_const, 2, mean)

R0_main_Cqui_infprob_const_peaks = sapply(apply(R0_main_Cqui_infprob_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_infprob_const_peaks_mean = mean(R0_main_Cqui_infprob_const_peaks)
R0_main_Cqui_infprob_const_peaks_0025 = quantile(R0_main_Cqui_infprob_const_peaks, probs=c(0.025))
R0_main_Cqui_infprob_const_peaks_0975 = quantile(R0_main_Cqui_infprob_const_peaks, probs=c(0.975))

R0_main_Cqui_infprob_const_Tmax = apply(R0_main_Cqui_infprob_const,1,FUN=Tmax)
R0_main_Cqui_infprob_const_Tmax_mean = mean(R0_main_Cqui_infprob_const_Tmax)
R0_main_Cqui_infprob_const_Tmax_0025 = quantile(R0_main_Cqui_infprob_const_Tmax, probs=c(0.025))
R0_main_Cqui_infprob_const_Tmax_0975 = quantile(R0_main_Cqui_infprob_const_Tmax, probs=c(0.975))

R0_main_Cqui_infprob_const_Tmin = apply(R0_main_Cqui_infprob_const,1,FUN=Tmin)
R0_main_Cqui_infprob_const_Tmin_mean = mean(R0_main_Cqui_infprob_const_Tmin)
R0_main_Cqui_infprob_const_Tmin_0025 = quantile(R0_main_Cqui_infprob_const_Tmin, probs=c(0.025))
R0_main_Cqui_infprob_const_Tmin_0975 = quantile(R0_main_Cqui_infprob_const_Tmin, probs=c(0.975))

# EIP

R0_main_Cqui_EIP_const = R0_main_f(biting_fit_Cqui, infprob_fit_pop, EIP_mean, 
                                   lf_fit_Cqui, omega, sur_fit_Cqui, ER, 
                                   egg_viability_fit_Cqui, dev_fit_Cqui)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cqui_EIP_const <= 0) == ncol(R0_main_Cqui_EIP_const))){
  print("only zeros")
  R0_main_Cqui_EIP_const <- R0_main_Cqui_EIP_const[-which(rowSums(R0_main_Cqui_EIP_const <= 0) == ncol(R0_main_Cqui_EIP_const)),]  
}

R0_main_Cqui_EIP_const_mean = apply(R0_main_Cqui_EIP_const, 2, mean)

R0_main_Cqui_EIP_const_peaks = sapply(apply(R0_main_Cqui_EIP_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_EIP_const_peaks_mean = mean(R0_main_Cqui_EIP_const_peaks)
R0_main_Cqui_EIP_const_peaks_0025 = quantile(R0_main_Cqui_EIP_const_peaks, probs=c(0.025))
R0_main_Cqui_EIP_const_peaks_0975 = quantile(R0_main_Cqui_EIP_const_peaks, probs=c(0.975))

R0_main_Cqui_EIP_const_Tmax = apply(R0_main_Cqui_EIP_const,1,FUN=Tmax)
R0_main_Cqui_EIP_const_Tmax_mean = mean(R0_main_Cqui_EIP_const_Tmax)
R0_main_Cqui_EIP_const_Tmax_0025 = quantile(R0_main_Cqui_EIP_const_Tmax, probs=c(0.025))
R0_main_Cqui_EIP_const_Tmax_0975 = quantile(R0_main_Cqui_EIP_const_Tmax, probs=c(0.975))

R0_main_Cqui_EIP_const_Tmin = apply(R0_main_Cqui_EIP_const,1,FUN=Tmin)
R0_main_Cqui_EIP_const_Tmin_mean = mean(R0_main_Cqui_EIP_const_Tmin)
R0_main_Cqui_EIP_const_Tmin_0025 = quantile(R0_main_Cqui_EIP_const_Tmin, probs=c(0.025))
R0_main_Cqui_EIP_const_Tmin_0975 = quantile(R0_main_Cqui_EIP_const_Tmin, probs=c(0.975))

# Adult lifespan

R0_main_Cqui_lf_const = R0_main_f(biting_fit_Cqui, infprob_fit_pop, EIP_fit_pop, 
                                  lf_data_mean, omega, sur_fit_Cqui, ER, 
                                  egg_viability_fit_Cqui, dev_fit_Cqui)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cqui_lf_const <= 0) == ncol(R0_main_Cqui_lf_const))){
  print("only zeros")
  R0_main_Cqui_lf_const <- R0_main_Cqui_lf_const[-which(rowSums(R0_main_Cqui_lf_const <= 0) == ncol(R0_main_Cqui_lf_const)),]  
}

R0_main_Cqui_lf_const_mean = apply(R0_main_Cqui_lf_const, 2, mean)

R0_main_Cqui_lf_const_peaks = sapply(apply(R0_main_Cqui_lf_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_lf_const_peaks_mean = mean(R0_main_Cqui_lf_const_peaks)
R0_main_Cqui_lf_const_peaks_0025 = quantile(R0_main_Cqui_lf_const_peaks, probs=c(0.025))
R0_main_Cqui_lf_const_peaks_0975 = quantile(R0_main_Cqui_lf_const_peaks, probs=c(0.975))

R0_main_Cqui_lf_const_Tmax = apply(R0_main_Cqui_lf_const,1,FUN=Tmax)
R0_main_Cqui_lf_const_Tmax_mean = mean(R0_main_Cqui_lf_const_Tmax)
R0_main_Cqui_lf_const_Tmax_0025 = quantile(R0_main_Cqui_lf_const_Tmax, probs=c(0.025))
R0_main_Cqui_lf_const_Tmax_0975 = quantile(R0_main_Cqui_lf_const_Tmax, probs=c(0.975))

R0_main_Cqui_lf_const_Tmin = apply(R0_main_Cqui_lf_const,1,FUN=Tmin)
R0_main_Cqui_lf_const_Tmin_mean = mean(R0_main_Cqui_lf_const_Tmin)
R0_main_Cqui_lf_const_Tmin_0025 = quantile(R0_main_Cqui_lf_const_Tmin, probs=c(0.025))
R0_main_Cqui_lf_const_Tmin_0975 = quantile(R0_main_Cqui_lf_const_Tmin, probs=c(0.975))

# Juvenile survival

R0_main_Cqui_sur_const = R0_main_f(biting_fit_Cqui, infprob_fit_pop, EIP_fit_pop, 
                                   lf_fit_Cqui, omega, sur_data_mean, ER, 
                                   egg_viability_fit_Cqui, dev_fit_Cqui)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cqui_sur_const <= 0) == ncol(R0_main_Cqui_sur_const))){
  print("only zeros")
  R0_main_Cqui_sur_const <- R0_main_Cqui_sur_const[-which(rowSums(R0_main_Cqui_sur_const <= 0) == ncol(R0_main_Cqui_sur_const)),]  
}

R0_main_Cqui_sur_const_mean = apply(R0_main_Cqui_sur_const, 2, mean)

R0_main_Cqui_sur_const_peaks = sapply(apply(R0_main_Cqui_sur_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_sur_const_peaks_mean = mean(R0_main_Cqui_sur_const_peaks)
R0_main_Cqui_sur_const_peaks_0025 = quantile(R0_main_Cqui_sur_const_peaks, probs=c(0.025))
R0_main_Cqui_sur_const_peaks_0975 = quantile(R0_main_Cqui_sur_const_peaks, probs=c(0.975))

R0_main_Cqui_sur_const_Tmax = apply(R0_main_Cqui_sur_const,1,FUN=Tmax)
R0_main_Cqui_sur_const_Tmax_mean = mean(R0_main_Cqui_sur_const_Tmax)
R0_main_Cqui_sur_const_Tmax_0025 = quantile(R0_main_Cqui_sur_const_Tmax, probs=c(0.025))
R0_main_Cqui_sur_const_Tmax_0975 = quantile(R0_main_Cqui_sur_const_Tmax, probs=c(0.975))

R0_main_Cqui_sur_const_Tmin = apply(R0_main_Cqui_sur_const,1,FUN=Tmin)
R0_main_Cqui_sur_const_Tmin_mean = mean(R0_main_Cqui_sur_const_Tmin)
R0_main_Cqui_sur_const_Tmin_0025 = quantile(R0_main_Cqui_sur_const_Tmin, probs=c(0.025))
R0_main_Cqui_sur_const_Tmin_0975 = quantile(R0_main_Cqui_sur_const_Tmin, probs=c(0.975))

# Egg viability

R0_main_Cqui_EV_const = R0_main_f(biting_fit_Cqui, infprob_fit_pop, EIP_fit_pop, 
                                  lf_fit_Cqui, omega, sur_fit_Cqui, ER, 
                                  egg_viability_data_mean, dev_fit_Cqui)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cqui_EV_const <= 0) == ncol(R0_main_Cqui_EV_const))){
  print("only zeros")
  R0_main_Cqui_EV_const <- R0_main_Cqui_EV_const[-which(rowSums(R0_main_Cqui_EV_const <= 0) == ncol(R0_main_Cqui_EV_const)),]  
}

R0_main_Cqui_EV_const_mean = apply(R0_main_Cqui_EV_const, 2, mean)

R0_main_Cqui_EV_const_peaks = sapply(apply(R0_main_Cqui_EV_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_EV_const_peaks_mean = mean(R0_main_Cqui_EV_const_peaks)
R0_main_Cqui_EV_const_peaks_0025 = quantile(R0_main_Cqui_EV_const_peaks, probs=c(0.025))
R0_main_Cqui_EV_const_peaks_0975 = quantile(R0_main_Cqui_EV_const_peaks, probs=c(0.975))

R0_main_Cqui_EV_const_Tmax = apply(R0_main_Cqui_EV_const,1,FUN=Tmax)
R0_main_Cqui_EV_const_Tmax_mean = mean(R0_main_Cqui_EV_const_Tmax)
R0_main_Cqui_EV_const_Tmax_0025 = quantile(R0_main_Cqui_EV_const_Tmax, probs=c(0.025))
R0_main_Cqui_EV_const_Tmax_0975 = quantile(R0_main_Cqui_EV_const_Tmax, probs=c(0.975))

R0_main_Cqui_EV_const_Tmin = apply(R0_main_Cqui_EV_const,1,FUN=Tmin)
R0_main_Cqui_EV_const_Tmin_mean = mean(R0_main_Cqui_EV_const_Tmin)
R0_main_Cqui_EV_const_Tmin_0025 = quantile(R0_main_Cqui_EV_const_Tmin, probs=c(0.025))
R0_main_Cqui_EV_const_Tmin_0975 = quantile(R0_main_Cqui_EV_const_Tmin, probs=c(0.975))

# Juvenile development rate

R0_main_Cqui_dev_const = R0_main_f(biting_fit_Cqui, infprob_fit_pop, EIP_fit_pop, 
                                   lf_fit_Cqui, omega, sur_fit_Cqui, ER, 
                                   egg_viability_fit_Cqui, dev_data_mean)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cqui_dev_const <= 0) == ncol(R0_main_Cqui_dev_const))){
  print("only zeros")
  R0_main_Cqui_dev_const <- R0_main_Cqui_dev_const[-which(rowSums(R0_main_Cqui_dev_const <= 0) == ncol(R0_main_Cqui_dev_const)),]  
}

R0_main_Cqui_dev_const_mean = apply(R0_main_Cqui_dev_const, 2, mean)

R0_main_Cqui_dev_const_peaks = sapply(apply(R0_main_Cqui_dev_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_dev_const_peaks_mean = mean(R0_main_Cqui_dev_const_peaks)
R0_main_Cqui_dev_const_peaks_0025 = quantile(R0_main_Cqui_dev_const_peaks, probs=c(0.025))
R0_main_Cqui_dev_const_peaks_0975 = quantile(R0_main_Cqui_dev_const_peaks, probs=c(0.975))

R0_main_Cqui_dev_const_Tmax = apply(R0_main_Cqui_dev_const,1,FUN=Tmax)
R0_main_Cqui_dev_const_Tmax_mean = mean(R0_main_Cqui_dev_const_Tmax)
R0_main_Cqui_dev_const_Tmax_0025 = quantile(R0_main_Cqui_dev_const_Tmax, probs=c(0.025))
R0_main_Cqui_dev_const_Tmax_0975 = quantile(R0_main_Cqui_dev_const_Tmax, probs=c(0.975))

R0_main_Cqui_dev_const_Tmin = apply(R0_main_Cqui_dev_const,1,FUN=Tmin)
R0_main_Cqui_dev_const_Tmin_mean = mean(R0_main_Cqui_dev_const_Tmin)
R0_main_Cqui_dev_const_Tmin_0025 = quantile(R0_main_Cqui_dev_const_Tmin, probs=c(0.025))
R0_main_Cqui_dev_const_Tmin_0975 = quantile(R0_main_Cqui_dev_const_Tmin, probs=c(0.975))

# Summarize results for Cx. quinquefasciatus in dataframe

df_R0_main_Cqui_const <- data.frame(x = temp,
                                    y = c(R0_main_Cqui_biting_const_mean/max(R0_main_Cqui_biting_const_mean),
                                          R0_main_Cqui_infprob_const_mean/max(R0_main_Cqui_infprob_const_mean),
                                          R0_main_Cqui_EIP_const_mean/max(R0_main_Cqui_EIP_const_mean),
                                          R0_main_Cqui_lf_const_mean/max(R0_main_Cqui_lf_const_mean),
                                          R0_main_Cqui_sur_const_mean/max(R0_main_Cqui_sur_const_mean),
                                          R0_main_Cqui_EV_const_mean/max(R0_main_Cqui_EV_const_mean),
                                          R0_main_Cqui_dev_const_mean/max(R0_main_Cqui_dev_const_mean),
                                          R0_main_Cqui_mean/max(R0_main_Cqui_mean)
                                    ),
                                    trait = c(rep("3Biting rate", length(temp)),
                                              rep("2Inf. prob.", length(temp)),
                                              rep("1EIP", length(temp)),
                                              rep("5Lifespan", length(temp)),
                                              rep("6Juv. survival", length(temp)),
                                              rep("4Egg viab.", length(temp)),
                                              rep("7Juv. dev. rate", length(temp)),
                                              rep("8None", length(temp))
                                    )
)

df_R0_main_Cqui_const_stats <- data.frame(peaks = c(R0_main_Cqui_biting_const_peaks_mean,
                                                    R0_main_Cqui_infprob_const_peaks_mean,
                                                    R0_main_Cqui_EIP_const_peaks_mean,
                                                    R0_main_Cqui_lf_const_peaks_mean,
                                                    R0_main_Cqui_sur_const_peaks_mean,
                                                    R0_main_Cqui_EV_const_peaks_mean,
                                                    R0_main_Cqui_dev_const_peaks_mean,
                                                    R0_main_Cqui_peaks_mean),
                                          peaks_lower = c(R0_main_Cqui_biting_const_peaks_0025,
                                                          R0_main_Cqui_infprob_const_peaks_0025,
                                                          R0_main_Cqui_EIP_const_peaks_0025,
                                                          R0_main_Cqui_lf_const_peaks_0025,
                                                          R0_main_Cqui_sur_const_peaks_0025,
                                                          R0_main_Cqui_EV_const_peaks_0025,
                                                          R0_main_Cqui_dev_const_peaks_0025,
                                                          R0_main_Cqui_peaks_0025),
                                          peaks_upper = c(R0_main_Cqui_biting_const_peaks_0975,
                                                          R0_main_Cqui_infprob_const_peaks_0975,
                                                          R0_main_Cqui_EIP_const_peaks_0975,
                                                          R0_main_Cqui_lf_const_peaks_0975,
                                                          R0_main_Cqui_sur_const_peaks_0975,
                                                          R0_main_Cqui_EV_const_peaks_0975,
                                                          R0_main_Cqui_dev_const_peaks_0975,
                                                          R0_main_Cqui_peaks_0975),
                                          Tmin = c(R0_main_Cqui_biting_const_Tmin_mean,
                                                   R0_main_Cqui_infprob_const_Tmin_mean,
                                                   R0_main_Cqui_EIP_const_Tmin_mean,
                                                   R0_main_Cqui_lf_const_Tmin_mean,
                                                   R0_main_Cqui_sur_const_Tmin_mean,
                                                   R0_main_Cqui_EV_const_Tmin_mean,
                                                   R0_main_Cqui_dev_const_Tmin_mean,
                                                   R0_main_Cqui_Tmin_mean),
                                          Tmin_lower = c(R0_main_Cqui_biting_const_Tmin_0025,
                                                         R0_main_Cqui_infprob_const_Tmin_0025,
                                                         R0_main_Cqui_EIP_const_Tmin_0025,
                                                         R0_main_Cqui_lf_const_Tmin_0025,
                                                         R0_main_Cqui_sur_const_Tmin_0025,
                                                         R0_main_Cqui_EV_const_Tmin_0025,
                                                         R0_main_Cqui_dev_const_Tmin_0025,
                                                         R0_main_Cqui_Tmin_0025),
                                          Tmin_upper = c(R0_main_Cqui_biting_const_Tmin_0975,
                                                         R0_main_Cqui_infprob_const_Tmin_0975,
                                                         R0_main_Cqui_EIP_const_Tmin_0975,
                                                         R0_main_Cqui_lf_const_Tmin_0975,
                                                         R0_main_Cqui_sur_const_Tmin_0975,
                                                         R0_main_Cqui_EV_const_Tmin_0975,
                                                         R0_main_Cqui_dev_const_Tmin_0975,
                                                         R0_main_Cqui_Tmin_0975),
                                          Tmax = c(R0_main_Cqui_biting_const_Tmax_mean,
                                                   R0_main_Cqui_infprob_const_Tmax_mean,
                                                   R0_main_Cqui_EIP_const_Tmax_mean,
                                                   R0_main_Cqui_lf_const_Tmax_mean,
                                                   R0_main_Cqui_sur_const_Tmax_mean,
                                                   R0_main_Cqui_EV_const_Tmax_mean,
                                                   R0_main_Cqui_dev_const_Tmax_mean,
                                                   R0_main_Cqui_Tmax_mean),
                                          Tmax_lower = c(R0_main_Cqui_biting_const_Tmax_0025,
                                                         R0_main_Cqui_infprob_const_Tmax_0025,
                                                         R0_main_Cqui_EIP_const_Tmax_0025,
                                                         R0_main_Cqui_lf_const_Tmax_0025,
                                                         R0_main_Cqui_sur_const_Tmax_0025,
                                                         R0_main_Cqui_EV_const_Tmax_0025,
                                                         R0_main_Cqui_dev_const_Tmax_0025,
                                                         R0_main_Cqui_Tmax_0025),
                                          Tmax_upper = c(R0_main_Cqui_biting_const_Tmax_0975,
                                                         R0_main_Cqui_infprob_const_Tmax_0975,
                                                         R0_main_Cqui_EIP_const_Tmax_0975,
                                                         R0_main_Cqui_lf_const_Tmax_0975,
                                                         R0_main_Cqui_sur_const_Tmax_0975,
                                                         R0_main_Cqui_EV_const_Tmax_0975,
                                                         R0_main_Cqui_dev_const_Tmax_0975,
                                                         R0_main_Cqui_Tmax_0975),
                                          trait = c("3Biting rate", "2Inf. prob.", "1EIP",
                                                    "5Lifespan", "6Juv. survival", "4Egg viab.", "7Juv. dev. rate",
                                                    "8None")
)

plot_Cqui_const <- ggplot() +
  geom_line(df_R0_main_Cqui_const, mapping = aes(x = x, y = y, color = trait), linewidth=0.6) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  theme_bw() +
  labs(x = "Temperature (°C)",
       color = "Constant trait", 
       title = expression("Mean temperature response of " * R[0]^rel)) +
  scale_color_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                                  "Lifespan", "Juv. survival", "Juv. dev. rate",
                                  "None")) +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_Cqui_stats <- ggplot(df_R0_main_Cqui_const_stats) +
  geom_pointrange(aes(x=peaks, xmin=peaks_lower, xmax=peaks_upper, y=trait, color=trait), linewidth=1.2, alpha=0.5, orientation = "y") +
  geom_pointrange(aes(x=Tmin, xmin=Tmin_lower, xmax=Tmin_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  geom_pointrange(aes(x=Tmax, xmin=Tmax_lower, xmax=Tmax_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  guides(color = "none") +
  theme_bw() +
  labs(x = "Temperature (°C)",
       title = "      Temperature limits and optimal temperature") +
  #scale_y_discrete(expand = c(0.1, 0.1)) +
  scale_y_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                              "Lifespan", "Juv. survival", "Juv. dev. rate",
                              "None")) +
  theme(plot.margin = unit(c(0.4, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cqui_const, plot_Cqui_stats, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. quinquefasciatus", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/R0_const_Cqui.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Sobol' analysis

r = 151:301
n = length(r)

sobols_first_Cqui <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

sobols_total_Cqui <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

for(i in 1:n){
  input_matrix <- data.frame(a = biting_fit_Cqui[,r[i]], bM = infprob_fit_pop[,r[i]], EIP = EIP_fit_pop[,r[i]], 
                             lf = lf_fit_Cqui[,r[i]], surJ = sur_fit_Cqui[,r[i]], 
                             EV= egg_viability_fit_Cqui[,r[i]], devJ = dev_fit_Cqui[,r[i]])
  
  X1 <- input_matrix[X_index1, ]
  X2 <- input_matrix[X_index2, ]
  
  sobol_result <- sobolmartinez(model = R0_main_sobol, 
                                X1 = X1, 
                                X2 = X2, 
                                nboot = 100)
  
  sobols_first_Cqui[i,-1] <- sobol_result$S$original
  sobols_total_Cqui[i,-1] <- sobol_result$T$original
}
sobols_first_Cqui[sobols_first_Cqui < 0] <- 0
sobols_total_Cqui[sobols_total_Cqui < 0] <- 0

sobols_first_long_Cqui <- sobols_first_Cqui %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cqui_sobol_first <- ggplot(sobols_first_long_Cqui) +
  geom_line(mapping = aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) + 
  labs(x = "Temperature (°C)",
       title = "First-order Sobol' indices",
       color = "Trait") +
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

sobols_total_long_Cqui <- sobols_total_Cqui %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cqui_sobol_total <-ggplot(sobols_total_long_Cqui) +
  geom_line(aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) +  
  labs(x = "Temperature (°C)",
       title = "Total-effect Sobol' indices",
       color = "Trait") + 
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cqui_sobol_first, plot_Cqui_sobol_total, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. quinquefasciatus", 
                                        size = 11, face = "bold.italic"))
plot

#ggsave("Figures/Sobol_Cqui.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Cx. pip. molestus

# calculate R0 samples from the trait samples
R0_main_Cmol = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                         lf_fit_Cmol, omega, sur_fit_Cmol, ER, 
                         egg_viability_fit_Cmol, dev_fit_Cmol)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cmol <= 0) == ncol(R0_main_Cmol))){
  print("only zeros")
  R0_main_Cmol <- R0_main_Cmol[-which(rowSums(R0_main_Cmol <= 0) == ncol(R0_main_Cmol)),]  
}

# calculate statistics of relative R0 at each temperature
R0_main_Cmol_mean = apply(R0_main_Cmol, 2, mean)
R0_main_Cmol_0025 = apply(R0_main_Cmol, 2, quantile, probs=c(0.025))
R0_main_Cmol_0975 = apply(R0_main_Cmol, 2, quantile, probs=c(0.975))

# calculate Topt for each sample and calculate statistics
R0_main_Cmol_peaks = sapply(apply(R0_main_Cmol,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_peaks_mean = mean(R0_main_Cmol_peaks)
R0_main_Cmol_peaks_median = median(R0_main_Cmol_peaks)
R0_main_Cmol_peaks_0025 = quantile(R0_main_Cmol_peaks, probs=c(0.025))
R0_main_Cmol_peaks_0975 = quantile(R0_main_Cmol_peaks, probs=c(0.975))

# calculate Tmin for each sample and calculate statistics
R0_main_Cmol_Tmin = apply(R0_main_Cmol,1,FUN=Tmin)
R0_main_Cmol_Tmin_mean = mean(R0_main_Cmol_Tmin)
R0_main_Cmol_Tmin_median = median(R0_main_Cmol_Tmin)
R0_main_Cmol_Tmin_0025 = quantile(R0_main_Cmol_Tmin, probs=c(0.025))
R0_main_Cmol_Tmin_0975 = quantile(R0_main_Cmol_Tmin, probs=c(0.975))

# calculate Tmax for each sample and calculate statistics
R0_main_Cmol_Tmax = apply(R0_main_Cmol,1,FUN=Tmax)
R0_main_Cmol_Tmax_mean = mean(R0_main_Cmol_Tmax)
R0_main_Cmol_Tmax_median = median(R0_main_Cmol_Tmax)
R0_main_Cmol_Tmax_0025 = quantile(R0_main_Cmol_Tmax, probs=c(0.025))
R0_main_Cmol_Tmax_0975 = quantile(R0_main_Cmol_Tmax, probs=c(0.975))

# calculate R0 samples again but leaving one trait constant

# Biting rate

R0_main_Cmol_biting_const = R0_main_f(biting_data_mean, infprob_fit_pop, EIP_fit_pop, 
                                      lf_fit_Cmol, omega, sur_fit_Cmol, ER, 
                                      egg_viability_fit_Cmol, dev_fit_Cmol)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cmol_biting_const <= 0) == ncol(R0_main_Cmol_biting_const))){
  print("only zeros")
  R0_main_Cmol_biting_const <- R0_main_Cmol_biting_const[-which(rowSums(R0_main_Cmol_biting_const <= 0) == ncol(R0_main_Cmol_biting_const)),]  
}

R0_main_Cmol_biting_const_mean = apply(R0_main_Cmol_biting_const, 2, mean)

R0_main_Cmol_biting_const_peaks = sapply(apply(R0_main_Cmol_biting_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_biting_const_peaks_mean = mean(R0_main_Cmol_biting_const_peaks)
R0_main_Cmol_biting_const_peaks_0025 = quantile(R0_main_Cmol_biting_const_peaks, probs=c(0.025))
R0_main_Cmol_biting_const_peaks_0975 = quantile(R0_main_Cmol_biting_const_peaks, probs=c(0.975))

R0_main_Cmol_biting_const_Tmax = apply(R0_main_Cmol_biting_const,1,FUN=Tmax)
R0_main_Cmol_biting_const_Tmax_mean = mean(R0_main_Cmol_biting_const_Tmax)
R0_main_Cmol_biting_const_Tmax_0025 = quantile(R0_main_Cmol_biting_const_Tmax, probs=c(0.025))
R0_main_Cmol_biting_const_Tmax_0975 = quantile(R0_main_Cmol_biting_const_Tmax, probs=c(0.975))

R0_main_Cmol_biting_const_Tmin = apply(R0_main_Cmol_biting_const,1,FUN=Tmin)
R0_main_Cmol_biting_const_Tmin_mean = mean(R0_main_Cmol_biting_const_Tmin)
R0_main_Cmol_biting_const_Tmin_0025 = quantile(R0_main_Cmol_biting_const_Tmin, probs=c(0.025))
R0_main_Cmol_biting_const_Tmin_0975 = quantile(R0_main_Cmol_biting_const_Tmin, probs=c(0.975))

# Infection probability

R0_main_Cmol_infprob_const = R0_main_f(biting_fit_pop, infprob_data_mean, EIP_fit_pop, 
                                       lf_fit_Cmol, omega, sur_fit_Cmol, ER, 
                                       egg_viability_fit_Cmol, dev_fit_Cmol)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cmol_infprob_const <= 0) == ncol(R0_main_Cmol_infprob_const))){
  print("only zeros")
  R0_main_Cmol_infprob_const <- R0_main_Cmol_infprob_const[-which(rowSums(R0_main_Cmol_infprob_const <= 0) == ncol(R0_main_Cmol_infprob_const)),]  
}

R0_main_Cmol_infprob_const_mean = apply(R0_main_Cmol_infprob_const, 2, mean)

R0_main_Cmol_infprob_const_peaks = sapply(apply(R0_main_Cmol_infprob_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_infprob_const_peaks_mean = mean(R0_main_Cmol_infprob_const_peaks)
R0_main_Cmol_infprob_const_peaks_0025 = quantile(R0_main_Cmol_infprob_const_peaks, probs=c(0.025))
R0_main_Cmol_infprob_const_peaks_0975 = quantile(R0_main_Cmol_infprob_const_peaks, probs=c(0.975))

R0_main_Cmol_infprob_const_Tmax = apply(R0_main_Cmol_infprob_const,1,FUN=Tmax)
R0_main_Cmol_infprob_const_Tmax_mean = mean(R0_main_Cmol_infprob_const_Tmax)
R0_main_Cmol_infprob_const_Tmax_0025 = quantile(R0_main_Cmol_infprob_const_Tmax, probs=c(0.025))
R0_main_Cmol_infprob_const_Tmax_0975 = quantile(R0_main_Cmol_infprob_const_Tmax, probs=c(0.975))

R0_main_Cmol_infprob_const_Tmin = apply(R0_main_Cmol_infprob_const,1,FUN=Tmin)
R0_main_Cmol_infprob_const_Tmin_mean = mean(R0_main_Cmol_infprob_const_Tmin)
R0_main_Cmol_infprob_const_Tmin_0025 = quantile(R0_main_Cmol_infprob_const_Tmin, probs=c(0.025))
R0_main_Cmol_infprob_const_Tmin_0975 = quantile(R0_main_Cmol_infprob_const_Tmin, probs=c(0.975))

# EIP

R0_main_Cmol_EIP_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_mean, 
                                   lf_fit_Cmol, omega, sur_fit_Cmol, ER, 
                                   egg_viability_fit_Cmol, dev_fit_Cmol)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cmol_EIP_const <= 0) == ncol(R0_main_Cmol_EIP_const))){
  print("only zeros")
  R0_main_Cmol_EIP_const <- R0_main_Cmol_EIP_const[-which(rowSums(R0_main_Cmol_EIP_const <= 0) == ncol(R0_main_Cmol_EIP_const)),]  
}

R0_main_Cmol_EIP_const_mean = apply(R0_main_Cmol_EIP_const, 2, mean)

R0_main_Cmol_EIP_const_peaks = sapply(apply(R0_main_Cmol_EIP_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_EIP_const_peaks_mean = mean(R0_main_Cmol_EIP_const_peaks)
R0_main_Cmol_EIP_const_peaks_0025 = quantile(R0_main_Cmol_EIP_const_peaks, probs=c(0.025))
R0_main_Cmol_EIP_const_peaks_0975 = quantile(R0_main_Cmol_EIP_const_peaks, probs=c(0.975))

R0_main_Cmol_EIP_const_Tmax = apply(R0_main_Cmol_EIP_const,1,FUN=Tmax)
R0_main_Cmol_EIP_const_Tmax_mean = mean(R0_main_Cmol_EIP_const_Tmax)
R0_main_Cmol_EIP_const_Tmax_0025 = quantile(R0_main_Cmol_EIP_const_Tmax, probs=c(0.025))
R0_main_Cmol_EIP_const_Tmax_0975 = quantile(R0_main_Cmol_EIP_const_Tmax, probs=c(0.975))

R0_main_Cmol_EIP_const_Tmin = apply(R0_main_Cmol_EIP_const,1,FUN=Tmin)
R0_main_Cmol_EIP_const_Tmin_mean = mean(R0_main_Cmol_EIP_const_Tmin)
R0_main_Cmol_EIP_const_Tmin_0025 = quantile(R0_main_Cmol_EIP_const_Tmin, probs=c(0.025))
R0_main_Cmol_EIP_const_Tmin_0975 = quantile(R0_main_Cmol_EIP_const_Tmin, probs=c(0.975))

# Adult lifespan

R0_main_Cmol_lf_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                                  lf_data_mean, omega, sur_fit_Cmol, ER, 
                                  egg_viability_fit_Cmol, dev_fit_Cmol)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cmol_lf_const <= 0) == ncol(R0_main_Cmol_lf_const))){
  print("only zeros")
  R0_main_Cmol_lf_const <- R0_main_Cmol_lf_const[-which(rowSums(R0_main_Cmol_lf_const <= 0) == ncol(R0_main_Cmol_lf_const)),]  
}

R0_main_Cmol_lf_const_mean = apply(R0_main_Cmol_lf_const, 2, mean)

R0_main_Cmol_lf_const_peaks = sapply(apply(R0_main_Cmol_lf_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_lf_const_peaks_mean = mean(R0_main_Cmol_lf_const_peaks)
R0_main_Cmol_lf_const_peaks_0025 = quantile(R0_main_Cmol_lf_const_peaks, probs=c(0.025))
R0_main_Cmol_lf_const_peaks_0975 = quantile(R0_main_Cmol_lf_const_peaks, probs=c(0.975))

R0_main_Cmol_lf_const_Tmax = apply(R0_main_Cmol_lf_const,1,FUN=Tmax)
R0_main_Cmol_lf_const_Tmax_mean = mean(R0_main_Cmol_lf_const_Tmax)
R0_main_Cmol_lf_const_Tmax_0025 = quantile(R0_main_Cmol_lf_const_Tmax, probs=c(0.025))
R0_main_Cmol_lf_const_Tmax_0975 = quantile(R0_main_Cmol_lf_const_Tmax, probs=c(0.975))

R0_main_Cmol_lf_const_Tmin = apply(R0_main_Cmol_lf_const,1,FUN=Tmin)
R0_main_Cmol_lf_const_Tmin_mean = mean(R0_main_Cmol_lf_const_Tmin)
R0_main_Cmol_lf_const_Tmin_0025 = quantile(R0_main_Cmol_lf_const_Tmin, probs=c(0.025))
R0_main_Cmol_lf_const_Tmin_0975 = quantile(R0_main_Cmol_lf_const_Tmin, probs=c(0.975))

# Juvenile survival

R0_main_Cmol_sur_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                                   lf_fit_Cmol, omega, sur_data_mean, ER, 
                                   egg_viability_fit_Cmol, dev_fit_Cmol)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cmol_sur_const <= 0) == ncol(R0_main_Cmol_sur_const))){
  print("only zeros")
  R0_main_Cmol_sur_const <- R0_main_Cmol_sur_const[-which(rowSums(R0_main_Cmol_sur_const <= 0) == ncol(R0_main_Cmol_sur_const)),]  
}

R0_main_Cmol_sur_const_mean = apply(R0_main_Cmol_sur_const, 2, mean)

R0_main_Cmol_sur_const_peaks = sapply(apply(R0_main_Cmol_sur_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_sur_const_peaks_mean = mean(R0_main_Cmol_sur_const_peaks)
R0_main_Cmol_sur_const_peaks_0025 = quantile(R0_main_Cmol_sur_const_peaks, probs=c(0.025))
R0_main_Cmol_sur_const_peaks_0975 = quantile(R0_main_Cmol_sur_const_peaks, probs=c(0.975))

R0_main_Cmol_sur_const_Tmax = apply(R0_main_Cmol_sur_const,1,FUN=Tmax)
R0_main_Cmol_sur_const_Tmax_mean = mean(R0_main_Cmol_sur_const_Tmax)
R0_main_Cmol_sur_const_Tmax_0025 = quantile(R0_main_Cmol_sur_const_Tmax, probs=c(0.025))
R0_main_Cmol_sur_const_Tmax_0975 = quantile(R0_main_Cmol_sur_const_Tmax, probs=c(0.975))

R0_main_Cmol_sur_const_Tmin = apply(R0_main_Cmol_sur_const,1,FUN=Tmin)
R0_main_Cmol_sur_const_Tmin_mean = mean(R0_main_Cmol_sur_const_Tmin)
R0_main_Cmol_sur_const_Tmin_0025 = quantile(R0_main_Cmol_sur_const_Tmin, probs=c(0.025))
R0_main_Cmol_sur_const_Tmin_0975 = quantile(R0_main_Cmol_sur_const_Tmin, probs=c(0.975))

# Egg viability

R0_main_Cmol_EV_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                                  lf_fit_Cmol, omega, sur_fit_Cmol, ER, 
                                  egg_viability_data_mean, dev_fit_Cmol)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cmol_EV_const <= 0) == ncol(R0_main_Cmol_EV_const))){
  print("only zeros")
  R0_main_Cmol_EV_const <- R0_main_Cmol_EV_const[-which(rowSums(R0_main_Cmol_EV_const <= 0) == ncol(R0_main_Cmol_EV_const)),]  
}

R0_main_Cmol_EV_const_mean = apply(R0_main_Cmol_EV_const, 2, mean)

R0_main_Cmol_EV_const_peaks = sapply(apply(R0_main_Cmol_EV_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_EV_const_peaks_mean = mean(R0_main_Cmol_EV_const_peaks)
R0_main_Cmol_EV_const_peaks_0025 = quantile(R0_main_Cmol_EV_const_peaks, probs=c(0.025))
R0_main_Cmol_EV_const_peaks_0975 = quantile(R0_main_Cmol_EV_const_peaks, probs=c(0.975))

R0_main_Cmol_EV_const_Tmax = apply(R0_main_Cmol_EV_const,1,FUN=Tmax)
R0_main_Cmol_EV_const_Tmax_mean = mean(R0_main_Cmol_EV_const_Tmax)
R0_main_Cmol_EV_const_Tmax_0025 = quantile(R0_main_Cmol_EV_const_Tmax, probs=c(0.025))
R0_main_Cmol_EV_const_Tmax_0975 = quantile(R0_main_Cmol_EV_const_Tmax, probs=c(0.975))

R0_main_Cmol_EV_const_Tmin = apply(R0_main_Cmol_EV_const,1,FUN=Tmin)
R0_main_Cmol_EV_const_Tmin_mean = mean(R0_main_Cmol_EV_const_Tmin)
R0_main_Cmol_EV_const_Tmin_0025 = quantile(R0_main_Cmol_EV_const_Tmin, probs=c(0.025))
R0_main_Cmol_EV_const_Tmin_0975 = quantile(R0_main_Cmol_EV_const_Tmin, probs=c(0.975))

# Juvenile development rate

R0_main_Cmol_dev_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                                   lf_fit_Cmol, omega, sur_fit_Cmol, ER, 
                                   egg_viability_fit_Cmol, dev_data_mean)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cmol_dev_const <= 0) == ncol(R0_main_Cmol_dev_const))){
  print("only zeros")
  R0_main_Cmol_dev_const <- R0_main_Cmol_dev_const[-which(rowSums(R0_main_Cmol_dev_const <= 0) == ncol(R0_main_Cmol_dev_const)),]  
}

R0_main_Cmol_dev_const_mean = apply(R0_main_Cmol_dev_const, 2, mean)

R0_main_Cmol_dev_const_peaks = sapply(apply(R0_main_Cmol_dev_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_dev_const_peaks_mean = mean(R0_main_Cmol_dev_const_peaks)
R0_main_Cmol_dev_const_peaks_0025 = quantile(R0_main_Cmol_dev_const_peaks, probs=c(0.025))
R0_main_Cmol_dev_const_peaks_0975 = quantile(R0_main_Cmol_dev_const_peaks, probs=c(0.975))

R0_main_Cmol_dev_const_Tmax = apply(R0_main_Cmol_dev_const,1,FUN=Tmax)
R0_main_Cmol_dev_const_Tmax_mean = mean(R0_main_Cmol_dev_const_Tmax)
R0_main_Cmol_dev_const_Tmax_0025 = quantile(R0_main_Cmol_dev_const_Tmax, probs=c(0.025))
R0_main_Cmol_dev_const_Tmax_0975 = quantile(R0_main_Cmol_dev_const_Tmax, probs=c(0.975))

R0_main_Cmol_dev_const_Tmin = apply(R0_main_Cmol_dev_const,1,FUN=Tmin)
R0_main_Cmol_dev_const_Tmin_mean = mean(R0_main_Cmol_dev_const_Tmin)
R0_main_Cmol_dev_const_Tmin_0025 = quantile(R0_main_Cmol_dev_const_Tmin, probs=c(0.025))
R0_main_Cmol_dev_const_Tmin_0975 = quantile(R0_main_Cmol_dev_const_Tmin, probs=c(0.975))

# Summarize results for Cx. pip. molestus in dataframe

df_R0_main_Cmol_const <- data.frame(x = temp,
                                    y = c(R0_main_Cmol_biting_const_mean/max(R0_main_Cmol_biting_const_mean),
                                          R0_main_Cmol_infprob_const_mean/max(R0_main_Cmol_infprob_const_mean),
                                          R0_main_Cmol_EIP_const_mean/max(R0_main_Cmol_EIP_const_mean),
                                          R0_main_Cmol_lf_const_mean/max(R0_main_Cmol_lf_const_mean),
                                          R0_main_Cmol_sur_const_mean/max(R0_main_Cmol_sur_const_mean),
                                          R0_main_Cmol_EV_const_mean/max(R0_main_Cmol_EV_const_mean),
                                          R0_main_Cmol_dev_const_mean/max(R0_main_Cmol_dev_const_mean),
                                          R0_main_Cmol_mean/max(R0_main_Cmol_mean)
                                    ),
                                    trait = c(rep("3Biting rate", length(temp)),
                                              rep("2Inf. prob.", length(temp)),
                                              rep("1EIP", length(temp)),
                                              rep("5Lifespan", length(temp)),
                                              rep("6Juv. survival", length(temp)),
                                              rep("4Egg viab.", length(temp)),
                                              rep("7Juv. dev. rate", length(temp)),
                                              rep("8None", length(temp))
                                    )
)

df_R0_main_Cmol_const_stats <- data.frame(peaks = c(R0_main_Cmol_biting_const_peaks_mean,
                                                    R0_main_Cmol_infprob_const_peaks_mean,
                                                    R0_main_Cmol_EIP_const_peaks_mean,
                                                    R0_main_Cmol_lf_const_peaks_mean,
                                                    R0_main_Cmol_sur_const_peaks_mean,
                                                    R0_main_Cmol_EV_const_peaks_mean,
                                                    R0_main_Cmol_dev_const_peaks_mean,
                                                    R0_main_Cmol_peaks_mean),
                                          peaks_lower = c(R0_main_Cmol_biting_const_peaks_0025,
                                                          R0_main_Cmol_infprob_const_peaks_0025,
                                                          R0_main_Cmol_EIP_const_peaks_0025,
                                                          R0_main_Cmol_lf_const_peaks_0025,
                                                          R0_main_Cmol_sur_const_peaks_0025,
                                                          R0_main_Cmol_EV_const_peaks_0025,
                                                          R0_main_Cmol_dev_const_peaks_0025,
                                                          R0_main_Cmol_peaks_0025),
                                          peaks_upper = c(R0_main_Cmol_biting_const_peaks_0975,
                                                          R0_main_Cmol_infprob_const_peaks_0975,
                                                          R0_main_Cmol_EIP_const_peaks_0975,
                                                          R0_main_Cmol_lf_const_peaks_0975,
                                                          R0_main_Cmol_sur_const_peaks_0975,
                                                          R0_main_Cmol_EV_const_peaks_0975,
                                                          R0_main_Cmol_dev_const_peaks_0975,
                                                          R0_main_Cmol_peaks_0975),
                                          Tmin = c(R0_main_Cmol_biting_const_Tmin_mean,
                                                   R0_main_Cmol_infprob_const_Tmin_mean,
                                                   R0_main_Cmol_EIP_const_Tmin_mean,
                                                   R0_main_Cmol_lf_const_Tmin_mean,
                                                   R0_main_Cmol_sur_const_Tmin_mean,
                                                   R0_main_Cmol_EV_const_Tmin_mean,
                                                   R0_main_Cmol_dev_const_Tmin_mean,
                                                   R0_main_Cmol_Tmin_mean),
                                          Tmin_lower = c(R0_main_Cmol_biting_const_Tmin_0025,
                                                         R0_main_Cmol_infprob_const_Tmin_0025,
                                                         R0_main_Cmol_EIP_const_Tmin_0025,
                                                         R0_main_Cmol_lf_const_Tmin_0025,
                                                         R0_main_Cmol_sur_const_Tmin_0025,
                                                         R0_main_Cmol_EV_const_Tmin_0025,
                                                         R0_main_Cmol_dev_const_Tmin_0025,
                                                         R0_main_Cmol_Tmin_0025),
                                          Tmin_upper = c(R0_main_Cmol_biting_const_Tmin_0975,
                                                         R0_main_Cmol_infprob_const_Tmin_0975,
                                                         R0_main_Cmol_EIP_const_Tmin_0975,
                                                         R0_main_Cmol_lf_const_Tmin_0975,
                                                         R0_main_Cmol_sur_const_Tmin_0975,
                                                         R0_main_Cmol_EV_const_Tmin_0975,
                                                         R0_main_Cmol_dev_const_Tmin_0975,
                                                         R0_main_Cmol_Tmin_0975),
                                          Tmax = c(R0_main_Cmol_biting_const_Tmax_mean,
                                                   R0_main_Cmol_infprob_const_Tmax_mean,
                                                   R0_main_Cmol_EIP_const_Tmax_mean,
                                                   R0_main_Cmol_lf_const_Tmax_mean,
                                                   R0_main_Cmol_sur_const_Tmax_mean,
                                                   R0_main_Cmol_EV_const_Tmax_mean,
                                                   R0_main_Cmol_dev_const_Tmax_mean,
                                                   R0_main_Cmol_Tmax_mean),
                                          Tmax_lower = c(R0_main_Cmol_biting_const_Tmax_0025,
                                                         R0_main_Cmol_infprob_const_Tmax_0025,
                                                         R0_main_Cmol_EIP_const_Tmax_0025,
                                                         R0_main_Cmol_lf_const_Tmax_0025,
                                                         R0_main_Cmol_sur_const_Tmax_0025,
                                                         R0_main_Cmol_EV_const_Tmax_0025,
                                                         R0_main_Cmol_dev_const_Tmax_0025,
                                                         R0_main_Cmol_Tmax_0025),
                                          Tmax_upper = c(R0_main_Cmol_biting_const_Tmax_0975,
                                                         R0_main_Cmol_infprob_const_Tmax_0975,
                                                         R0_main_Cmol_EIP_const_Tmax_0975,
                                                         R0_main_Cmol_lf_const_Tmax_0975,
                                                         R0_main_Cmol_sur_const_Tmax_0975,
                                                         R0_main_Cmol_EV_const_Tmax_0975,
                                                         R0_main_Cmol_dev_const_Tmax_0975,
                                                         R0_main_Cmol_Tmax_0975),
                                          trait = c("3Biting rate", "2Inf. prob.", "1EIP",
                                                    "5Lifespan", "6Juv. survival", "4Egg viab.", "7Juv. dev. rate",
                                                    "8None")
)

plot_Cmol_const <- ggplot() +
  geom_line(df_R0_main_Cmol_const, mapping = aes(x = x, y = y, color = trait), linewidth=0.6) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. pipiens molestus")))) +
  labs(x = "Temperature (°C)",
       color = "Constant trait", 
       title = expression("Mean temperature response of " * R[0]^rel)) +
  scale_color_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                                  "Lifespan", "Juv. survival", "Juv. dev. rate",
                                  "None")) +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_Cmol_stats <- ggplot(df_R0_main_Cmol_const_stats) +
  geom_pointrange(aes(x=peaks, xmin=peaks_lower, xmax=peaks_upper, y=trait, color=trait), linewidth=1.2, alpha=0.5, orientation = "y") +
  geom_pointrange(aes(x=Tmin, xmin=Tmin_lower, xmax=Tmin_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  geom_pointrange(aes(x=Tmax, xmin=Tmax_lower, xmax=Tmax_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  guides(color = "none") +
  theme_bw() +
  labs(x = "Temperature (°C)",
       title = "      Temperature limits and optimal temperature") +
  #scale_y_discrete(expand = c(0.1, 0.1)) +
  scale_y_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                              "Lifespan", "Juv. survival", "Juv. dev. rate",
                              "None")) +
  theme(plot.margin = unit(c(0.4, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cmol_const, plot_Cmol_stats, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. pipiens molestus", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/R0_const_Cmol.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Sobol' analysis

r = 151:301
n = length(r)

sobols_first_Cmol <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

sobols_total_Cmol <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

for(i in 1:n){
  input_matrix <- data.frame(a = biting_fit_pop[,r[i]], bM = infprob_fit_pop[,r[i]], EIP = EIP_fit_pop[,r[i]], 
                             lf = lf_fit_Cmol[,r[i]], surJ = sur_fit_Cmol[,r[i]], 
                             EV = egg_viability_fit_Cmol[,r[i]], devJ = dev_fit_Cmol[,r[i]])
  
  X1 <- input_matrix[X_index1, ]
  X2 <- input_matrix[X_index2, ]
  
  sobol_result <- sobolmartinez(model = R0_main_sobol, 
                                X1 = X1, 
                                X2 = X2, 
                                nboot = 100)
  
  sobols_first_Cmol[i,-1] <- sobol_result$S$original
  sobols_total_Cmol[i,-1] <- sobol_result$T$original
}
sobols_first_Cmol[sobols_first_Cmol < 0] <- 0
sobols_total_Cmol[sobols_total_Cmol < 0] <- 0

sobols_first_long_Cmol <- sobols_first_Cmol %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cmol_sobol_first <- ggplot(sobols_first_long_Cmol) +
  geom_line(mapping = aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) + 
  labs(x = "Temperature (°C)",
       title = "First-order Sobol' indices",
       color = "Trait") +
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

sobols_total_long_Cmol <- sobols_total_Cmol %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cmol_sobol_total <-ggplot(sobols_total_long_Cmol) +
  geom_line(aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) +  
  labs(x = "Temperature (°C)",
       title = "Total-effect Sobol' indices",
       color = "Trait") + 
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cmol_sobol_first, plot_Cmol_sobol_total, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. pipiens molestus", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/Sobol_Cmol.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Cx. pip. pallens

# calculate R0 samples from the trait samples
R0_main_Cpal = R0_main_f(biting_fit_Cpal, infprob_fit_pop, EIP_fit_pop, 
                         lf_fit_Cpal, omega, sur_fit_Cpal, ER, 
                         egg_viability_fit_Cpal, dev_fit_Cpal)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal <= 0) == ncol(R0_main_Cpal))){
  print("only zeros")
  R0_main_Cpal <- R0_main_Cpal[-which(rowSums(R0_main_Cpal <= 0) == ncol(R0_main_Cpal)),]  
}

# calculate statistics of relative R0 at each temperature
R0_main_Cpal_mean = apply(R0_main_Cpal, 2, mean)
R0_main_Cpal_0025 = apply(R0_main_Cpal, 2, quantile, probs=c(0.025))
R0_main_Cpal_0975 = apply(R0_main_Cpal, 2, quantile, probs=c(0.975))

# calculate Topt for each sample and calculate statistics
R0_main_Cpal_peaks = sapply(apply(R0_main_Cpal,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_peaks_mean = mean(R0_main_Cpal_peaks)
R0_main_Cpal_peaks_median = median(R0_main_Cpal_peaks)
R0_main_Cpal_peaks_0025 = quantile(R0_main_Cpal_peaks, probs=c(0.025))
R0_main_Cpal_peaks_0975 = quantile(R0_main_Cpal_peaks, probs=c(0.975))

# calculate Tmin for each sample and calculate statistics
R0_main_Cpal_Tmin = apply(R0_main_Cpal,1,FUN=Tmin)
R0_main_Cpal_Tmin_mean = mean(R0_main_Cpal_Tmin)
R0_main_Cpal_Tmin_median = median(R0_main_Cpal_Tmin)
R0_main_Cpal_Tmin_0025 = quantile(R0_main_Cpal_Tmin, probs=c(0.025))
R0_main_Cpal_Tmin_0975 = quantile(R0_main_Cpal_Tmin, probs=c(0.975))

# calculate Tmax for each sample and calculate statistics
R0_main_Cpal_Tmax = apply(R0_main_Cpal,1,FUN=Tmax)
R0_main_Cpal_Tmax_mean = mean(R0_main_Cpal_Tmax)
R0_main_Cpal_Tmax_median = median(R0_main_Cpal_Tmax)
R0_main_Cpal_Tmax_0025 = quantile(R0_main_Cpal_Tmax, probs=c(0.025))
R0_main_Cpal_Tmax_0975 = quantile(R0_main_Cpal_Tmax, probs=c(0.975))

# calculate R0 samples again but leaving one trait constant

# Biting rate

R0_main_Cpal_biting_const = R0_main_f(biting_data_mean, infprob_fit_pop, EIP_fit_pop, 
                                      lf_fit_Cpal, omega, sur_fit_Cpal, ER, 
                                      egg_viability_fit_Cpal, dev_fit_Cpal)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal_biting_const <= 0) == ncol(R0_main_Cpal_biting_const))){
  print("only zeros")
  R0_main_Cpal_biting_const <- R0_main_Cpal_biting_const[-which(rowSums(R0_main_Cpal_biting_const <= 0) == ncol(R0_main_Cpal_biting_const)),]  
}

R0_main_Cpal_biting_const_mean = apply(R0_main_Cpal_biting_const, 2, mean)

R0_main_Cpal_biting_const_peaks = sapply(apply(R0_main_Cpal_biting_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_biting_const_peaks_mean = mean(R0_main_Cpal_biting_const_peaks)
R0_main_Cpal_biting_const_peaks_0025 = quantile(R0_main_Cpal_biting_const_peaks, probs=c(0.025))
R0_main_Cpal_biting_const_peaks_0975 = quantile(R0_main_Cpal_biting_const_peaks, probs=c(0.975))

R0_main_Cpal_biting_const_Tmax = apply(R0_main_Cpal_biting_const,1,FUN=Tmax)
R0_main_Cpal_biting_const_Tmax_mean = mean(R0_main_Cpal_biting_const_Tmax)
R0_main_Cpal_biting_const_Tmax_0025 = quantile(R0_main_Cpal_biting_const_Tmax, probs=c(0.025))
R0_main_Cpal_biting_const_Tmax_0975 = quantile(R0_main_Cpal_biting_const_Tmax, probs=c(0.975))

R0_main_Cpal_biting_const_Tmin = apply(R0_main_Cpal_biting_const,1,FUN=Tmin)
R0_main_Cpal_biting_const_Tmin_mean = mean(R0_main_Cpal_biting_const_Tmin)
R0_main_Cpal_biting_const_Tmin_0025 = quantile(R0_main_Cpal_biting_const_Tmin, probs=c(0.025))
R0_main_Cpal_biting_const_Tmin_0975 = quantile(R0_main_Cpal_biting_const_Tmin, probs=c(0.975))

# Infection probability

R0_main_Cpal_infprob_const = R0_main_f(biting_fit_Cpal, infprob_data_mean, EIP_fit_pop, 
                                       lf_fit_Cpal, omega, sur_fit_Cpal, ER, 
                                       egg_viability_fit_Cpal, dev_fit_Cpal)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal_infprob_const <= 0) == ncol(R0_main_Cpal_infprob_const))){
  print("only zeros")
  R0_main_Cpal_infprob_const <- R0_main_Cpal_infprob_const[-which(rowSums(R0_main_Cpal_infprob_const <= 0) == ncol(R0_main_Cpal_infprob_const)),]  
}

R0_main_Cpal_infprob_const_mean = apply(R0_main_Cpal_infprob_const, 2, mean)

R0_main_Cpal_infprob_const_peaks = sapply(apply(R0_main_Cpal_infprob_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_infprob_const_peaks_mean = mean(R0_main_Cpal_infprob_const_peaks)
R0_main_Cpal_infprob_const_peaks_0025 = quantile(R0_main_Cpal_infprob_const_peaks, probs=c(0.025))
R0_main_Cpal_infprob_const_peaks_0975 = quantile(R0_main_Cpal_infprob_const_peaks, probs=c(0.975))

R0_main_Cpal_infprob_const_Tmax = apply(R0_main_Cpal_infprob_const,1,FUN=Tmax)
R0_main_Cpal_infprob_const_Tmax_mean = mean(R0_main_Cpal_infprob_const_Tmax)
R0_main_Cpal_infprob_const_Tmax_0025 = quantile(R0_main_Cpal_infprob_const_Tmax, probs=c(0.025))
R0_main_Cpal_infprob_const_Tmax_0975 = quantile(R0_main_Cpal_infprob_const_Tmax, probs=c(0.975))

R0_main_Cpal_infprob_const_Tmin = apply(R0_main_Cpal_infprob_const,1,FUN=Tmin)
R0_main_Cpal_infprob_const_Tmin_mean = mean(R0_main_Cpal_infprob_const_Tmin)
R0_main_Cpal_infprob_const_Tmin_0025 = quantile(R0_main_Cpal_infprob_const_Tmin, probs=c(0.025))
R0_main_Cpal_infprob_const_Tmin_0975 = quantile(R0_main_Cpal_infprob_const_Tmin, probs=c(0.975))

# EIP

R0_main_Cpal_EIP_const = R0_main_f(biting_fit_Cpal, infprob_fit_pop, EIP_mean, 
                                   lf_fit_Cpal, omega, sur_fit_Cpal, ER, 
                                   egg_viability_fit_Cpal, dev_fit_Cpal)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal_EIP_const <= 0) == ncol(R0_main_Cpal_EIP_const))){
  print("only zeros")
  R0_main_Cpal_EIP_const <- R0_main_Cpal_EIP_const[-which(rowSums(R0_main_Cpal_EIP_const <= 0) == ncol(R0_main_Cpal_EIP_const)),]  
}

R0_main_Cpal_EIP_const_mean = apply(R0_main_Cpal_EIP_const, 2, mean)

R0_main_Cpal_EIP_const_peaks = sapply(apply(R0_main_Cpal_EIP_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_EIP_const_peaks_mean = mean(R0_main_Cpal_EIP_const_peaks)
R0_main_Cpal_EIP_const_peaks_0025 = quantile(R0_main_Cpal_EIP_const_peaks, probs=c(0.025))
R0_main_Cpal_EIP_const_peaks_0975 = quantile(R0_main_Cpal_EIP_const_peaks, probs=c(0.975))

R0_main_Cpal_EIP_const_Tmax = apply(R0_main_Cpal_EIP_const,1,FUN=Tmax)
R0_main_Cpal_EIP_const_Tmax_mean = mean(R0_main_Cpal_EIP_const_Tmax)
R0_main_Cpal_EIP_const_Tmax_0025 = quantile(R0_main_Cpal_EIP_const_Tmax, probs=c(0.025))
R0_main_Cpal_EIP_const_Tmax_0975 = quantile(R0_main_Cpal_EIP_const_Tmax, probs=c(0.975))

R0_main_Cpal_EIP_const_Tmin = apply(R0_main_Cpal_EIP_const,1,FUN=Tmin)
R0_main_Cpal_EIP_const_Tmin_mean = mean(R0_main_Cpal_EIP_const_Tmin)
R0_main_Cpal_EIP_const_Tmin_0025 = quantile(R0_main_Cpal_EIP_const_Tmin, probs=c(0.025))
R0_main_Cpal_EIP_const_Tmin_0975 = quantile(R0_main_Cpal_EIP_const_Tmin, probs=c(0.975))

# Adult lifespan

R0_main_Cpal_lf_const = R0_main_f(biting_fit_Cpal, infprob_fit_pop, EIP_fit_pop, 
                                  lf_data_mean, omega, sur_fit_Cpal, ER, 
                                  egg_viability_fit_Cpal, dev_fit_Cpal)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal_lf_const <= 0) == ncol(R0_main_Cpal_lf_const))){
  print("only zeros")
  R0_main_Cpal_lf_const <- R0_main_Cpal_lf_const[-which(rowSums(R0_main_Cpal_lf_const <= 0) == ncol(R0_main_Cpal_lf_const)),]  
}

R0_main_Cpal_lf_const_mean = apply(R0_main_Cpal_lf_const, 2, mean)

R0_main_Cpal_lf_const_peaks = sapply(apply(R0_main_Cpal_lf_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_lf_const_peaks_mean = mean(R0_main_Cpal_lf_const_peaks)
R0_main_Cpal_lf_const_peaks_0025 = quantile(R0_main_Cpal_lf_const_peaks, probs=c(0.025))
R0_main_Cpal_lf_const_peaks_0975 = quantile(R0_main_Cpal_lf_const_peaks, probs=c(0.975))

R0_main_Cpal_lf_const_Tmax = apply(R0_main_Cpal_lf_const,1,FUN=Tmax)
R0_main_Cpal_lf_const_Tmax_mean = mean(R0_main_Cpal_lf_const_Tmax)
R0_main_Cpal_lf_const_Tmax_0025 = quantile(R0_main_Cpal_lf_const_Tmax, probs=c(0.025))
R0_main_Cpal_lf_const_Tmax_0975 = quantile(R0_main_Cpal_lf_const_Tmax, probs=c(0.975))

R0_main_Cpal_lf_const_Tmin = apply(R0_main_Cpal_lf_const,1,FUN=Tmin)
R0_main_Cpal_lf_const_Tmin_mean = mean(R0_main_Cpal_lf_const_Tmin)
R0_main_Cpal_lf_const_Tmin_0025 = quantile(R0_main_Cpal_lf_const_Tmin, probs=c(0.025))
R0_main_Cpal_lf_const_Tmin_0975 = quantile(R0_main_Cpal_lf_const_Tmin, probs=c(0.975))

# Juvenile survival

R0_main_Cpal_sur_const = R0_main_f(biting_fit_Cpal, infprob_fit_pop, EIP_fit_pop, 
                                   lf_fit_Cpal, omega, sur_data_mean, ER, 
                                   egg_viability_fit_Cpal, dev_fit_Cpal)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal_sur_const <= 0) == ncol(R0_main_Cpal_sur_const))){
  print("only zeros")
  R0_main_Cpal_sur_const <- R0_main_Cpal_sur_const[-which(rowSums(R0_main_Cpal_sur_const <= 0) == ncol(R0_main_Cpal_sur_const)),]  
}

R0_main_Cpal_sur_const_mean = apply(R0_main_Cpal_sur_const, 2, mean)

R0_main_Cpal_sur_const_peaks = sapply(apply(R0_main_Cpal_sur_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_sur_const_peaks_mean = mean(R0_main_Cpal_sur_const_peaks)
R0_main_Cpal_sur_const_peaks_0025 = quantile(R0_main_Cpal_sur_const_peaks, probs=c(0.025))
R0_main_Cpal_sur_const_peaks_0975 = quantile(R0_main_Cpal_sur_const_peaks, probs=c(0.975))

R0_main_Cpal_sur_const_Tmax = apply(R0_main_Cpal_sur_const,1,FUN=Tmax)
R0_main_Cpal_sur_const_Tmax_mean = mean(R0_main_Cpal_sur_const_Tmax)
R0_main_Cpal_sur_const_Tmax_0025 = quantile(R0_main_Cpal_sur_const_Tmax, probs=c(0.025))
R0_main_Cpal_sur_const_Tmax_0975 = quantile(R0_main_Cpal_sur_const_Tmax, probs=c(0.975))

R0_main_Cpal_sur_const_Tmin = apply(R0_main_Cpal_sur_const,1,FUN=Tmin)
R0_main_Cpal_sur_const_Tmin_mean = mean(R0_main_Cpal_sur_const_Tmin)
R0_main_Cpal_sur_const_Tmin_0025 = quantile(R0_main_Cpal_sur_const_Tmin, probs=c(0.025))
R0_main_Cpal_sur_const_Tmin_0975 = quantile(R0_main_Cpal_sur_const_Tmin, probs=c(0.975))

# Egg viability

R0_main_Cpal_EV_const = R0_main_f(biting_fit_Cpal, infprob_fit_pop, EIP_fit_pop, 
                                  lf_fit_Cpal, omega, sur_fit_Cpal, ER, 
                                  egg_viability_data_mean, dev_fit_Cpal)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal_EV_const <= 0) == ncol(R0_main_Cpal_EV_const))){
  print("only zeros")
  R0_main_Cpal_EV_const <- R0_main_Cpal_EV_const[-which(rowSums(R0_main_Cpal_EV_const <= 0) == ncol(R0_main_Cpal_EV_const)),]  
}

R0_main_Cpal_EV_const_mean = apply(R0_main_Cpal_EV_const, 2, mean)

R0_main_Cpal_EV_const_peaks = sapply(apply(R0_main_Cpal_EV_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_EV_const_peaks_mean = mean(R0_main_Cpal_EV_const_peaks)
R0_main_Cpal_EV_const_peaks_0025 = quantile(R0_main_Cpal_EV_const_peaks, probs=c(0.025))
R0_main_Cpal_EV_const_peaks_0975 = quantile(R0_main_Cpal_EV_const_peaks, probs=c(0.975))

R0_main_Cpal_EV_const_Tmax = apply(R0_main_Cpal_EV_const,1,FUN=Tmax)
R0_main_Cpal_EV_const_Tmax_mean = mean(R0_main_Cpal_EV_const_Tmax)
R0_main_Cpal_EV_const_Tmax_0025 = quantile(R0_main_Cpal_EV_const_Tmax, probs=c(0.025))
R0_main_Cpal_EV_const_Tmax_0975 = quantile(R0_main_Cpal_EV_const_Tmax, probs=c(0.975))

R0_main_Cpal_EV_const_Tmin = apply(R0_main_Cpal_EV_const,1,FUN=Tmin)
R0_main_Cpal_EV_const_Tmin_mean = mean(R0_main_Cpal_EV_const_Tmin)
R0_main_Cpal_EV_const_Tmin_0025 = quantile(R0_main_Cpal_EV_const_Tmin, probs=c(0.025))
R0_main_Cpal_EV_const_Tmin_0975 = quantile(R0_main_Cpal_EV_const_Tmin, probs=c(0.975))

# Juvenile development rate

R0_main_Cpal_dev_const = R0_main_f(biting_fit_Cpal, infprob_fit_pop, EIP_fit_pop, 
                                   lf_fit_Cpal, omega, sur_fit_Cpal, ER, 
                                   egg_viability_fit_Cpal, dev_data_mean)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal_dev_const <= 0) == ncol(R0_main_Cpal_dev_const))){
  print("only zeros")
  R0_main_Cpal_dev_const <- R0_main_Cpal_dev_const[-which(rowSums(R0_main_Cpal_dev_const <= 0) == ncol(R0_main_Cpal_dev_const)),]  
}

R0_main_Cpal_dev_const_mean = apply(R0_main_Cpal_dev_const, 2, mean)

R0_main_Cpal_dev_const_peaks = sapply(apply(R0_main_Cpal_dev_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_dev_const_peaks_mean = mean(R0_main_Cpal_dev_const_peaks)
R0_main_Cpal_dev_const_peaks_0025 = quantile(R0_main_Cpal_dev_const_peaks, probs=c(0.025))
R0_main_Cpal_dev_const_peaks_0975 = quantile(R0_main_Cpal_dev_const_peaks, probs=c(0.975))

R0_main_Cpal_dev_const_Tmax = apply(R0_main_Cpal_dev_const,1,FUN=Tmax)
R0_main_Cpal_dev_const_Tmax_mean = mean(R0_main_Cpal_dev_const_Tmax)
R0_main_Cpal_dev_const_Tmax_0025 = quantile(R0_main_Cpal_dev_const_Tmax, probs=c(0.025))
R0_main_Cpal_dev_const_Tmax_0975 = quantile(R0_main_Cpal_dev_const_Tmax, probs=c(0.975))

R0_main_Cpal_dev_const_Tmin = apply(R0_main_Cpal_dev_const,1,FUN=Tmin)
R0_main_Cpal_dev_const_Tmin_mean = mean(R0_main_Cpal_dev_const_Tmin)
R0_main_Cpal_dev_const_Tmin_0025 = quantile(R0_main_Cpal_dev_const_Tmin, probs=c(0.025))
R0_main_Cpal_dev_const_Tmin_0975 = quantile(R0_main_Cpal_dev_const_Tmin, probs=c(0.975))

# Summarize results for Cx. pip. molestus in dataframe

df_R0_main_Cpal_const <- data.frame(x = temp,
                                    y = c(R0_main_Cpal_biting_const_mean/max(R0_main_Cpal_biting_const_mean),
                                          R0_main_Cpal_infprob_const_mean/max(R0_main_Cpal_infprob_const_mean),
                                          R0_main_Cpal_EIP_const_mean/max(R0_main_Cpal_EIP_const_mean),
                                          R0_main_Cpal_lf_const_mean/max(R0_main_Cpal_lf_const_mean),
                                          R0_main_Cpal_sur_const_mean/max(R0_main_Cpal_sur_const_mean),
                                          R0_main_Cpal_EV_const_mean/max(R0_main_Cpal_EV_const_mean),
                                          R0_main_Cpal_dev_const_mean/max(R0_main_Cpal_dev_const_mean),
                                          R0_main_Cpal_mean/max(R0_main_Cpal_mean)
                                    ),
                                    trait = c(rep("3Biting rate", length(temp)),
                                              rep("2Inf. prob.", length(temp)),
                                              rep("1EIP", length(temp)),
                                              rep("5Lifespan", length(temp)),
                                              rep("6Juv. survival", length(temp)),
                                              rep("4Egg viab.", length(temp)),
                                              rep("7Juv. dev. rate", length(temp)),
                                              rep("8None", length(temp))
                                    )
)

df_R0_main_Cpal_const_stats <- data.frame(peaks = c(R0_main_Cpal_biting_const_peaks_mean,
                                                    R0_main_Cpal_infprob_const_peaks_mean,
                                                    R0_main_Cpal_EIP_const_peaks_mean,
                                                    R0_main_Cpal_lf_const_peaks_mean,
                                                    R0_main_Cpal_sur_const_peaks_mean,
                                                    R0_main_Cpal_EV_const_peaks_mean,
                                                    R0_main_Cpal_dev_const_peaks_mean,
                                                    R0_main_Cpal_peaks_mean),
                                          peaks_lower = c(R0_main_Cpal_biting_const_peaks_0025,
                                                          R0_main_Cpal_infprob_const_peaks_0025,
                                                          R0_main_Cpal_EIP_const_peaks_0025,
                                                          R0_main_Cpal_lf_const_peaks_0025,
                                                          R0_main_Cpal_sur_const_peaks_0025,
                                                          R0_main_Cpal_EV_const_peaks_0025,
                                                          R0_main_Cpal_dev_const_peaks_0025,
                                                          R0_main_Cpal_peaks_0025),
                                          peaks_upper = c(R0_main_Cpal_biting_const_peaks_0975,
                                                          R0_main_Cpal_infprob_const_peaks_0975,
                                                          R0_main_Cpal_EIP_const_peaks_0975,
                                                          R0_main_Cpal_lf_const_peaks_0975,
                                                          R0_main_Cpal_sur_const_peaks_0975,
                                                          R0_main_Cpal_EV_const_peaks_0975,
                                                          R0_main_Cpal_dev_const_peaks_0975,
                                                          R0_main_Cpal_peaks_0975),
                                          Tmin = c(R0_main_Cpal_biting_const_Tmin_mean,
                                                   R0_main_Cpal_infprob_const_Tmin_mean,
                                                   R0_main_Cpal_EIP_const_Tmin_mean,
                                                   R0_main_Cpal_lf_const_Tmin_mean,
                                                   R0_main_Cpal_sur_const_Tmin_mean,
                                                   R0_main_Cpal_EV_const_Tmin_mean,
                                                   R0_main_Cpal_dev_const_Tmin_mean,
                                                   R0_main_Cpal_Tmin_mean),
                                          Tmin_lower = c(R0_main_Cpal_biting_const_Tmin_0025,
                                                         R0_main_Cpal_infprob_const_Tmin_0025,
                                                         R0_main_Cpal_EIP_const_Tmin_0025,
                                                         R0_main_Cpal_lf_const_Tmin_0025,
                                                         R0_main_Cpal_sur_const_Tmin_0025,
                                                         R0_main_Cpal_EV_const_Tmin_0025,
                                                         R0_main_Cpal_dev_const_Tmin_0025,
                                                         R0_main_Cpal_Tmin_0025),
                                          Tmin_upper = c(R0_main_Cpal_biting_const_Tmin_0975,
                                                         R0_main_Cpal_infprob_const_Tmin_0975,
                                                         R0_main_Cpal_EIP_const_Tmin_0975,
                                                         R0_main_Cpal_lf_const_Tmin_0975,
                                                         R0_main_Cpal_sur_const_Tmin_0975,
                                                         R0_main_Cpal_EV_const_Tmin_0975,
                                                         R0_main_Cpal_dev_const_Tmin_0975,
                                                         R0_main_Cpal_Tmin_0975),
                                          Tmax = c(R0_main_Cpal_biting_const_Tmax_mean,
                                                   R0_main_Cpal_infprob_const_Tmax_mean,
                                                   R0_main_Cpal_EIP_const_Tmax_mean,
                                                   R0_main_Cpal_lf_const_Tmax_mean,
                                                   R0_main_Cpal_sur_const_Tmax_mean,
                                                   R0_main_Cpal_EV_const_Tmax_mean,
                                                   R0_main_Cpal_dev_const_Tmax_mean,
                                                   R0_main_Cpal_Tmax_mean),
                                          Tmax_lower = c(R0_main_Cpal_biting_const_Tmax_0025,
                                                         R0_main_Cpal_infprob_const_Tmax_0025,
                                                         R0_main_Cpal_EIP_const_Tmax_0025,
                                                         R0_main_Cpal_lf_const_Tmax_0025,
                                                         R0_main_Cpal_sur_const_Tmax_0025,
                                                         R0_main_Cpal_EV_const_Tmax_0025,
                                                         R0_main_Cpal_dev_const_Tmax_0025,
                                                         R0_main_Cpal_Tmax_0025),
                                          Tmax_upper = c(R0_main_Cpal_biting_const_Tmax_0975,
                                                         R0_main_Cpal_infprob_const_Tmax_0975,
                                                         R0_main_Cpal_EIP_const_Tmax_0975,
                                                         R0_main_Cpal_lf_const_Tmax_0975,
                                                         R0_main_Cpal_sur_const_Tmax_0975,
                                                         R0_main_Cpal_EV_const_Tmax_0975,
                                                         R0_main_Cpal_dev_const_Tmax_0975,
                                                         R0_main_Cpal_Tmax_0975),
                                          trait = c("3Biting rate", "2Inf. prob.", "1EIP",
                                                    "5Lifespan", "6Juv. survival", "4Egg viab.", "7Juv. dev. rate",
                                                    "8None")
)

plot_Cpal_const <- ggplot() +
  geom_line(df_R0_main_Cpal_const, mapping = aes(x = x, y = y, color = trait), linewidth=0.6) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. pipiens molestus")))) +
  labs(x = "Temperature (°C)",
       color = "Constant trait", 
       title = expression("Mean temperature response of " * R[0]^rel)) +
  scale_color_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                                  "Lifespan", "Juv. survival", "Juv. dev. rate",
                                  "None")) +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_Cpal_stats <- ggplot(df_R0_main_Cpal_const_stats) +
  geom_pointrange(aes(x=peaks, xmin=peaks_lower, xmax=peaks_upper, y=trait, color=trait), linewidth=1.2, alpha=0.5, orientation = "y") +
  geom_pointrange(aes(x=Tmin, xmin=Tmin_lower, xmax=Tmin_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  geom_pointrange(aes(x=Tmax, xmin=Tmax_lower, xmax=Tmax_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  guides(color = "none") +
  theme_bw() +
  labs(x = "Temperature (°C)",
       title = "      Temperature limits and optimal temperature") +
  scale_y_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                              "Lifespan", "Juv. survival", "Juv. dev. rate",
                              "None")) +
  theme(plot.margin = unit(c(0.4, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cpal_const, plot_Cpal_stats, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. pipiens pallens", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/R0_const_Cpal.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Sobol' analysis

r = 151:301
n = length(r)

sobols_first_Cpal <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

sobols_total_Cpal <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

for(i in 1:n){
  input_matrix <- data.frame(a = biting_fit_Cpal[,r[i]], bM = infprob_fit_pop[,r[i]], EIP = EIP_fit_pop[,r[i]], 
                             lf = lf_fit_Cpal[,r[i]], surJ = sur_fit_Cpal[,r[i]], 
                             EV = egg_viability_fit_Cpal[,r[i]], devJ = dev_fit_Cpal[,r[i]])
  
  X1 <- input_matrix[X_index1, ]
  X2 <- input_matrix[X_index2, ]
  
  sobol_result <- sobolmartinez(model = R0_main_sobol, 
                                X1 = X1, 
                                X2 = X2, 
                                nboot = 100)
  
  sobols_first_Cpal[i,-1] <- sobol_result$S$original
  sobols_total_Cpal[i,-1] <- sobol_result$T$original
}
sobols_first_Cpal[sobols_first_Cpal < 0] <- 0
sobols_total_Cpal[sobols_total_Cpal < 0] <- 0

sobols_first_long_Cpal <- sobols_first_Cpal %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cpal_sobol_first <- ggplot(sobols_first_long_Cpal) +
  geom_line(mapping = aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) + 
  labs(x = "Temperature (°C)",
       title = "First-order Sobol' indices",
       color = "Trait") +
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

sobols_total_long_Cpal <- sobols_total_Cpal %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cpal_sobol_total <-ggplot(sobols_total_long_Cpal) +
  geom_line(aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) +  
  labs(x = "Temperature (°C)",
       title = "Total-effect Sobol' indices",
       color = "Trait") + 
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cpal_sobol_first, plot_Cpal_sobol_total, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. pipiens pallens", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/Sobol_Cpal.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Cx. restuans

# calculate R0 samples from the trait samples
R0_main_Cres = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                         lf_fit_Cres, omega, sur_fit_Cres, ER, 
                         egg_viability_fit_pop, dev_fit_Cres)


# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres <= 0) == ncol(R0_main_Cres))){
  print("only zeros")
  R0_main_Cres <- R0_main_Cres[-which(rowSums(R0_main_Cres <= 0) == ncol(R0_main_Cres)),]  
}

# calculate statistics of relative R0 at each temperature
R0_main_Cres_mean = apply(R0_main_Cres, 2, mean)
R0_main_Cres_0025 = apply(R0_main_Cres, 2, quantile, probs=c(0.025))
R0_main_Cres_0975 = apply(R0_main_Cres, 2, quantile, probs=c(0.975))

# calculate Topt for each sample and calculate statistics
R0_main_Cres_peaks = sapply(apply(R0_main_Cres,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_peaks_mean = mean(R0_main_Cres_peaks)
R0_main_Cres_peaks_median = median(R0_main_Cres_peaks)
R0_main_Cres_peaks_0025 = quantile(R0_main_Cres_peaks, probs=c(0.025))
R0_main_Cres_peaks_0975 = quantile(R0_main_Cres_peaks, probs=c(0.975))

# calculate Tmin for each sample and calculate statistics
R0_main_Cres_Tmin = apply(R0_main_Cres,1,FUN=Tmin)
R0_main_Cres_Tmin_mean = mean(R0_main_Cres_Tmin)
R0_main_Cres_Tmin_median = median(R0_main_Cres_Tmin)
R0_main_Cres_Tmin_0025 = quantile(R0_main_Cres_Tmin, probs=c(0.025))
R0_main_Cres_Tmin_0975 = quantile(R0_main_Cres_Tmin, probs=c(0.975))

# calculate Tmax for each sample and calculate statistics
R0_main_Cres_Tmax = apply(R0_main_Cres,1,FUN=Tmax)
R0_main_Cres_Tmax_mean = mean(R0_main_Cres_Tmax)
R0_main_Cres_Tmax_median = median(R0_main_Cres_Tmax)
R0_main_Cres_Tmax_0025 = quantile(R0_main_Cres_Tmax, probs=c(0.025))
R0_main_Cres_Tmax_0975 = quantile(R0_main_Cres_Tmax, probs=c(0.975))

# calculate R0 samples again but leaving one trait constant

# Biting rate

R0_main_Cres_biting_const = R0_main_f(biting_data_mean, infprob_fit_pop, EIP_fit_pop, 
                                      lf_fit_Cres, omega, sur_fit_Cres, ER, 
                                      egg_viability_fit_pop, dev_fit_Cres)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres_biting_const <= 0) == ncol(R0_main_Cres_biting_const))){
  print("only zeros")
  R0_main_Cres_biting_const <- R0_main_Cres_biting_const[-which(rowSums(R0_main_Cres_biting_const <= 0) == ncol(R0_main_Cres_biting_const)),]  
}

R0_main_Cres_biting_const_mean = apply(R0_main_Cres_biting_const, 2, mean)

R0_main_Cres_biting_const_peaks = sapply(apply(R0_main_Cres_biting_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_biting_const_peaks_mean = mean(R0_main_Cres_biting_const_peaks)
R0_main_Cres_biting_const_peaks_0025 = quantile(R0_main_Cres_biting_const_peaks, probs=c(0.025))
R0_main_Cres_biting_const_peaks_0975 = quantile(R0_main_Cres_biting_const_peaks, probs=c(0.975))

R0_main_Cres_biting_const_Tmax = apply(R0_main_Cres_biting_const,1,FUN=Tmax)
R0_main_Cres_biting_const_Tmax_mean = mean(R0_main_Cres_biting_const_Tmax)
R0_main_Cres_biting_const_Tmax_0025 = quantile(R0_main_Cres_biting_const_Tmax, probs=c(0.025))
R0_main_Cres_biting_const_Tmax_0975 = quantile(R0_main_Cres_biting_const_Tmax, probs=c(0.975))

R0_main_Cres_biting_const_Tmin = apply(R0_main_Cres_biting_const,1,FUN=Tmin)
R0_main_Cres_biting_const_Tmin_mean = mean(R0_main_Cres_biting_const_Tmin)
R0_main_Cres_biting_const_Tmin_0025 = quantile(R0_main_Cres_biting_const_Tmin, probs=c(0.025))
R0_main_Cres_biting_const_Tmin_0975 = quantile(R0_main_Cres_biting_const_Tmin, probs=c(0.975))

# Infection probability

R0_main_Cres_infprob_const = R0_main_f(biting_fit_pop, infprob_data_mean, EIP_fit_pop, 
                                       lf_fit_Cres, omega, sur_fit_Cres, ER, 
                                       egg_viability_fit_pop, dev_fit_Cres)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres_infprob_const <= 0) == ncol(R0_main_Cres_infprob_const))){
  print("only zeros")
  R0_main_Cres_infprob_const <- R0_main_Cres_infprob_const[-which(rowSums(R0_main_Cres_infprob_const <= 0) == ncol(R0_main_Cres_infprob_const)),]  
}

R0_main_Cres_infprob_const_mean = apply(R0_main_Cres_infprob_const, 2, mean)

R0_main_Cres_infprob_const_peaks = sapply(apply(R0_main_Cres_infprob_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_infprob_const_peaks_mean = mean(R0_main_Cres_infprob_const_peaks)
R0_main_Cres_infprob_const_peaks_0025 = quantile(R0_main_Cres_infprob_const_peaks, probs=c(0.025))
R0_main_Cres_infprob_const_peaks_0975 = quantile(R0_main_Cres_infprob_const_peaks, probs=c(0.975))

R0_main_Cres_infprob_const_Tmax = apply(R0_main_Cres_infprob_const,1,FUN=Tmax)
R0_main_Cres_infprob_const_Tmax_mean = mean(R0_main_Cres_infprob_const_Tmax)
R0_main_Cres_infprob_const_Tmax_0025 = quantile(R0_main_Cres_infprob_const_Tmax, probs=c(0.025))
R0_main_Cres_infprob_const_Tmax_0975 = quantile(R0_main_Cres_infprob_const_Tmax, probs=c(0.975))

R0_main_Cres_infprob_const_Tmin = apply(R0_main_Cres_infprob_const,1,FUN=Tmin)
R0_main_Cres_infprob_const_Tmin_mean = mean(R0_main_Cres_infprob_const_Tmin)
R0_main_Cres_infprob_const_Tmin_0025 = quantile(R0_main_Cres_infprob_const_Tmin, probs=c(0.025))
R0_main_Cres_infprob_const_Tmin_0975 = quantile(R0_main_Cres_infprob_const_Tmin, probs=c(0.975))

# EIP

R0_main_Cres_EIP_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_mean, 
                                   lf_fit_Cres, omega, sur_fit_Cres, ER, 
                                   egg_viability_fit_pop, dev_fit_Cres)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres_EIP_const <= 0) == ncol(R0_main_Cres_EIP_const))){
  print("only zeros")
  R0_main_Cres_EIP_const <- R0_main_Cres_EIP_const[-which(rowSums(R0_main_Cres_EIP_const <= 0) == ncol(R0_main_Cres_EIP_const)),]  
}

R0_main_Cres_EIP_const_mean = apply(R0_main_Cres_EIP_const, 2, mean)

R0_main_Cres_EIP_const_peaks = sapply(apply(R0_main_Cres_EIP_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_EIP_const_peaks_mean = mean(R0_main_Cres_EIP_const_peaks)
R0_main_Cres_EIP_const_peaks_0025 = quantile(R0_main_Cres_EIP_const_peaks, probs=c(0.025))
R0_main_Cres_EIP_const_peaks_0975 = quantile(R0_main_Cres_EIP_const_peaks, probs=c(0.975))

R0_main_Cres_EIP_const_Tmax = apply(R0_main_Cres_EIP_const,1,FUN=Tmax)
R0_main_Cres_EIP_const_Tmax_mean = mean(R0_main_Cres_EIP_const_Tmax)
R0_main_Cres_EIP_const_Tmax_0025 = quantile(R0_main_Cres_EIP_const_Tmax, probs=c(0.025))
R0_main_Cres_EIP_const_Tmax_0975 = quantile(R0_main_Cres_EIP_const_Tmax, probs=c(0.975))

R0_main_Cres_EIP_const_Tmin = apply(R0_main_Cres_EIP_const,1,FUN=Tmin)
R0_main_Cres_EIP_const_Tmin_mean = mean(R0_main_Cres_EIP_const_Tmin)
R0_main_Cres_EIP_const_Tmin_0025 = quantile(R0_main_Cres_EIP_const_Tmin, probs=c(0.025))
R0_main_Cres_EIP_const_Tmin_0975 = quantile(R0_main_Cres_EIP_const_Tmin, probs=c(0.975))

# Adult lifespan

R0_main_Cres_lf_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                                  lf_data_mean, omega, sur_fit_Cres, ER, 
                                  egg_viability_fit_pop, dev_fit_Cres)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres_lf_const <= 0) == ncol(R0_main_Cres_lf_const))){
  print("only zeros")
  R0_main_Cres_lf_const <- R0_main_Cres_lf_const[-which(rowSums(R0_main_Cres_lf_const <= 0) == ncol(R0_main_Cres_lf_const)),]  
}

R0_main_Cres_lf_const_mean = apply(R0_main_Cres_lf_const, 2, mean)

R0_main_Cres_lf_const_peaks = sapply(apply(R0_main_Cres_lf_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_lf_const_peaks_mean = mean(R0_main_Cres_lf_const_peaks)
R0_main_Cres_lf_const_peaks_0025 = quantile(R0_main_Cres_lf_const_peaks, probs=c(0.025))
R0_main_Cres_lf_const_peaks_0975 = quantile(R0_main_Cres_lf_const_peaks, probs=c(0.975))

R0_main_Cres_lf_const_Tmax = apply(R0_main_Cres_lf_const,1,FUN=Tmax)
R0_main_Cres_lf_const_Tmax_mean = mean(R0_main_Cres_lf_const_Tmax)
R0_main_Cres_lf_const_Tmax_0025 = quantile(R0_main_Cres_lf_const_Tmax, probs=c(0.025))
R0_main_Cres_lf_const_Tmax_0975 = quantile(R0_main_Cres_lf_const_Tmax, probs=c(0.975))

R0_main_Cres_lf_const_Tmin = apply(R0_main_Cres_lf_const,1,FUN=Tmin)
R0_main_Cres_lf_const_Tmin_mean = mean(R0_main_Cres_lf_const_Tmin)
R0_main_Cres_lf_const_Tmin_0025 = quantile(R0_main_Cres_lf_const_Tmin, probs=c(0.025))
R0_main_Cres_lf_const_Tmin_0975 = quantile(R0_main_Cres_lf_const_Tmin, probs=c(0.975))

# Juvenile survival

R0_main_Cres_sur_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                                   lf_fit_Cres, omega, sur_data_mean, ER, 
                                   egg_viability_fit_pop, dev_fit_Cres)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres_sur_const <= 0) == ncol(R0_main_Cres_sur_const))){
  print("only zeros")
  R0_main_Cres_sur_const <- R0_main_Cres_sur_const[-which(rowSums(R0_main_Cres_sur_const <= 0) == ncol(R0_main_Cres_sur_const)),]  
}

R0_main_Cres_sur_const_mean = apply(R0_main_Cres_sur_const, 2, mean)

R0_main_Cres_sur_const_peaks = sapply(apply(R0_main_Cres_sur_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_sur_const_peaks_mean = mean(R0_main_Cres_sur_const_peaks)
R0_main_Cres_sur_const_peaks_0025 = quantile(R0_main_Cres_sur_const_peaks, probs=c(0.025))
R0_main_Cres_sur_const_peaks_0975 = quantile(R0_main_Cres_sur_const_peaks, probs=c(0.975))

R0_main_Cres_sur_const_Tmax = apply(R0_main_Cres_sur_const,1,FUN=Tmax)
R0_main_Cres_sur_const_Tmax_mean = mean(R0_main_Cres_sur_const_Tmax)
R0_main_Cres_sur_const_Tmax_0025 = quantile(R0_main_Cres_sur_const_Tmax, probs=c(0.025))
R0_main_Cres_sur_const_Tmax_0975 = quantile(R0_main_Cres_sur_const_Tmax, probs=c(0.975))

R0_main_Cres_sur_const_Tmin = apply(R0_main_Cres_sur_const,1,FUN=Tmin)
R0_main_Cres_sur_const_Tmin_mean = mean(R0_main_Cres_sur_const_Tmin)
R0_main_Cres_sur_const_Tmin_0025 = quantile(R0_main_Cres_sur_const_Tmin, probs=c(0.025))
R0_main_Cres_sur_const_Tmin_0975 = quantile(R0_main_Cres_sur_const_Tmin, probs=c(0.975))

# Egg viability

R0_main_Cres_EV_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                                  lf_fit_Cres, omega, sur_fit_Cres, ER, 
                                  egg_viability_data_mean, dev_fit_Cres)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres_EV_const <= 0) == ncol(R0_main_Cres_EV_const))){
  print("only zeros")
  R0_main_Cres_EV_const <- R0_main_Cres_EV_const[-which(rowSums(R0_main_Cres_EV_const <= 0) == ncol(R0_main_Cres_EV_const)),]  
}

R0_main_Cres_EV_const_mean = apply(R0_main_Cres_EV_const, 2, mean)

R0_main_Cres_EV_const_peaks = sapply(apply(R0_main_Cres_EV_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_EV_const_peaks_mean = mean(R0_main_Cres_EV_const_peaks)
R0_main_Cres_EV_const_peaks_0025 = quantile(R0_main_Cres_EV_const_peaks, probs=c(0.025))
R0_main_Cres_EV_const_peaks_0975 = quantile(R0_main_Cres_EV_const_peaks, probs=c(0.975))

R0_main_Cres_EV_const_Tmax = apply(R0_main_Cres_EV_const,1,FUN=Tmax)
R0_main_Cres_EV_const_Tmax_mean = mean(R0_main_Cres_EV_const_Tmax)
R0_main_Cres_EV_const_Tmax_0025 = quantile(R0_main_Cres_EV_const_Tmax, probs=c(0.025))
R0_main_Cres_EV_const_Tmax_0975 = quantile(R0_main_Cres_EV_const_Tmax, probs=c(0.975))

R0_main_Cres_EV_const_Tmin = apply(R0_main_Cres_EV_const,1,FUN=Tmin)
R0_main_Cres_EV_const_Tmin_mean = mean(R0_main_Cres_EV_const_Tmin)
R0_main_Cres_EV_const_Tmin_0025 = quantile(R0_main_Cres_EV_const_Tmin, probs=c(0.025))
R0_main_Cres_EV_const_Tmin_0975 = quantile(R0_main_Cres_EV_const_Tmin, probs=c(0.975))

# Juvenile development

R0_main_Cres_dev_const = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, 
                                   lf_fit_Cres, omega, sur_fit_Cres, ER, 
                                   egg_viability_fit_pop, dev_data_mean)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres_dev_const <= 0) == ncol(R0_main_Cres_dev_const))){
  print("only zeros")
  R0_main_Cres_dev_const <- R0_main_Cres_dev_const[-which(rowSums(R0_main_Cres_dev_const <= 0) == ncol(R0_main_Cres_dev_const)),]  
}

R0_main_Cres_dev_const_mean = apply(R0_main_Cres_dev_const, 2, mean)

R0_main_Cres_dev_const_peaks = sapply(apply(R0_main_Cres_dev_const,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_dev_const_peaks_mean = mean(R0_main_Cres_dev_const_peaks)
R0_main_Cres_dev_const_peaks_0025 = quantile(R0_main_Cres_dev_const_peaks, probs=c(0.025))
R0_main_Cres_dev_const_peaks_0975 = quantile(R0_main_Cres_dev_const_peaks, probs=c(0.975))

R0_main_Cres_dev_const_Tmax = apply(R0_main_Cres_dev_const,1,FUN=Tmax)
R0_main_Cres_dev_const_Tmax_mean = mean(R0_main_Cres_dev_const_Tmax)
R0_main_Cres_dev_const_Tmax_0025 = quantile(R0_main_Cres_dev_const_Tmax, probs=c(0.025))
R0_main_Cres_dev_const_Tmax_0975 = quantile(R0_main_Cres_dev_const_Tmax, probs=c(0.975))

R0_main_Cres_dev_const_Tmin = apply(R0_main_Cres_dev_const,1,FUN=Tmin)
R0_main_Cres_dev_const_Tmin_mean = mean(R0_main_Cres_dev_const_Tmin)
R0_main_Cres_dev_const_Tmin_0025 = quantile(R0_main_Cres_dev_const_Tmin, probs=c(0.025))
R0_main_Cres_dev_const_Tmin_0975 = quantile(R0_main_Cres_dev_const_Tmin, probs=c(0.975))

# Summarize results for Cx. pip. molestus in dataframe

df_R0_main_Cres_const <- data.frame(x = temp,
                                    y = c(R0_main_Cres_biting_const_mean/max(R0_main_Cres_biting_const_mean),
                                          R0_main_Cres_infprob_const_mean/max(R0_main_Cres_infprob_const_mean),
                                          R0_main_Cres_EIP_const_mean/max(R0_main_Cres_EIP_const_mean),
                                          R0_main_Cres_lf_const_mean/max(R0_main_Cres_lf_const_mean),
                                          R0_main_Cres_sur_const_mean/max(R0_main_Cres_sur_const_mean),
                                          R0_main_Cres_EV_const_mean/max(R0_main_Cres_EV_const_mean),
                                          R0_main_Cres_dev_const_mean/max(R0_main_Cres_dev_const_mean),
                                          R0_main_Cres_mean/max(R0_main_Cres_mean)
                                    ),
                                    trait = c(rep("3Biting rate", length(temp)),
                                              rep("2Inf. prob.", length(temp)),
                                              rep("1EIP", length(temp)),
                                              rep("5Lifespan", length(temp)),
                                              rep("6Juv. survival", length(temp)),
                                              rep("4Egg viab.", length(temp)),
                                              rep("7Juv. dev. rate", length(temp)),
                                              rep("8None", length(temp))
                                    )
)

df_R0_main_Cres_const_stats <- data.frame(peaks = c(R0_main_Cres_biting_const_peaks_mean,
                                                    R0_main_Cres_infprob_const_peaks_mean,
                                                    R0_main_Cres_EIP_const_peaks_mean,
                                                    R0_main_Cres_lf_const_peaks_mean,
                                                    R0_main_Cres_sur_const_peaks_mean,
                                                    R0_main_Cres_EV_const_peaks_mean,
                                                    R0_main_Cres_dev_const_peaks_mean,
                                                    R0_main_Cres_peaks_mean),
                                          peaks_lower = c(R0_main_Cres_biting_const_peaks_0025,
                                                          R0_main_Cres_infprob_const_peaks_0025,
                                                          R0_main_Cres_EIP_const_peaks_0025,
                                                          R0_main_Cres_lf_const_peaks_0025,
                                                          R0_main_Cres_sur_const_peaks_0025,
                                                          R0_main_Cres_EV_const_peaks_0025,
                                                          R0_main_Cres_dev_const_peaks_0025,
                                                          R0_main_Cres_peaks_0025),
                                          peaks_upper = c(R0_main_Cres_biting_const_peaks_0975,
                                                          R0_main_Cres_infprob_const_peaks_0975,
                                                          R0_main_Cres_EIP_const_peaks_0975,
                                                          R0_main_Cres_lf_const_peaks_0975,
                                                          R0_main_Cres_sur_const_peaks_0975,
                                                          R0_main_Cres_EV_const_peaks_0975,
                                                          R0_main_Cres_dev_const_peaks_0975,
                                                          R0_main_Cres_peaks_0975),
                                          Tmin = c(R0_main_Cres_biting_const_Tmin_mean,
                                                   R0_main_Cres_infprob_const_Tmin_mean,
                                                   R0_main_Cres_EIP_const_Tmin_mean,
                                                   R0_main_Cres_lf_const_Tmin_mean,
                                                   R0_main_Cres_sur_const_Tmin_mean,
                                                   R0_main_Cres_EV_const_Tmin_mean,
                                                   R0_main_Cres_dev_const_Tmin_mean,
                                                   R0_main_Cres_Tmin_mean),
                                          Tmin_lower = c(R0_main_Cres_biting_const_Tmin_0025,
                                                         R0_main_Cres_infprob_const_Tmin_0025,
                                                         R0_main_Cres_EIP_const_Tmin_0025,
                                                         R0_main_Cres_lf_const_Tmin_0025,
                                                         R0_main_Cres_sur_const_Tmin_0025,
                                                         R0_main_Cres_EV_const_Tmin_0025,
                                                         R0_main_Cres_dev_const_Tmin_0025,
                                                         R0_main_Cres_Tmin_0025),
                                          Tmin_upper = c(R0_main_Cres_biting_const_Tmin_0975,
                                                         R0_main_Cres_infprob_const_Tmin_0975,
                                                         R0_main_Cres_EIP_const_Tmin_0975,
                                                         R0_main_Cres_lf_const_Tmin_0975,
                                                         R0_main_Cres_sur_const_Tmin_0975,
                                                         R0_main_Cres_EV_const_Tmin_0975,
                                                         R0_main_Cres_dev_const_Tmin_0975,
                                                         R0_main_Cres_Tmin_0975),
                                          Tmax = c(R0_main_Cres_biting_const_Tmax_mean,
                                                   R0_main_Cres_infprob_const_Tmax_mean,
                                                   R0_main_Cres_EIP_const_Tmax_mean,
                                                   R0_main_Cres_lf_const_Tmax_mean,
                                                   R0_main_Cres_sur_const_Tmax_mean,
                                                   R0_main_Cres_EV_const_Tmax_mean,
                                                   R0_main_Cres_dev_const_Tmax_mean,
                                                   R0_main_Cres_Tmax_mean),
                                          Tmax_lower = c(R0_main_Cres_biting_const_Tmax_0025,
                                                         R0_main_Cres_infprob_const_Tmax_0025,
                                                         R0_main_Cres_EIP_const_Tmax_0025,
                                                         R0_main_Cres_lf_const_Tmax_0025,
                                                         R0_main_Cres_sur_const_Tmax_0025,
                                                         R0_main_Cres_EV_const_Tmax_0025,
                                                         R0_main_Cres_dev_const_Tmax_0025,
                                                         R0_main_Cres_Tmax_0025),
                                          Tmax_upper = c(R0_main_Cres_biting_const_Tmax_0975,
                                                         R0_main_Cres_infprob_const_Tmax_0975,
                                                         R0_main_Cres_EIP_const_Tmax_0975,
                                                         R0_main_Cres_lf_const_Tmax_0975,
                                                         R0_main_Cres_sur_const_Tmax_0975,
                                                         R0_main_Cres_EV_const_Tmax_0975,
                                                         R0_main_Cres_dev_const_Tmax_0975,
                                                         R0_main_Cres_Tmax_0975),
                                          trait = c("3Biting rate", "2Inf. prob.", "1EIP",
                                                    "5Lifespan", "6Juv. survival", "4Egg viab.", "7Juv. dev. rate",
                                                    "8None")
)

plot_Cres_const <- ggplot() +
  geom_line(df_R0_main_Cres_const, mapping = aes(x = x, y = y, color = trait), linewidth=0.6) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. restuans")))) +
  labs(x = "Temperature (°C)",
       color = "Constant trait", 
       title = expression("Mean temperature response of " * R[0]^rel)) +
  scale_color_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                                  "Lifespan", "Juv. survival", "Juv. dev. rate",
                                  "None")) +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_Cres_stats <- ggplot(df_R0_main_Cres_const_stats) +
  geom_pointrange(aes(x=peaks, xmin=peaks_lower, xmax=peaks_upper, y=trait, color=trait), linewidth=1.2, alpha=0.5, orientation = "y") +
  geom_pointrange(aes(x=Tmin, xmin=Tmin_lower, xmax=Tmin_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  geom_pointrange(aes(x=Tmax, xmin=Tmax_lower, xmax=Tmax_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  guides(color = "none") +
  theme_bw() +
  labs(x = "Temperature (°C)",
       title = "      Temperature limits and optimal temperature") +
  scale_y_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                              "Lifespan", "Juv. survival", "Juv. dev. rate",
                              "None")) +
  theme(plot.margin = unit(c(0.4, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cres_const, plot_Cres_stats, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. restuans", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/R0_const_Cres.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Sobol' analysis

r = 151:301
n = length(r)

sobols_first_Cres <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

sobols_total_Cres <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

for(i in 1:n){
  input_matrix <- data.frame(a = biting_fit_pop[,r[i]], bM = infprob_fit_pop[,r[i]], EIP = EIP_fit_pop[,r[i]], 
                             lf = lf_fit_Cres[,r[i]], surJ = sur_fit_Cres[,r[i]], 
                             EV= egg_viability_fit_pop[,r[i]], devJ = dev_fit_Cres[,r[i]])
  
  X1 <- input_matrix[X_index1, ]
  X2 <- input_matrix[X_index2, ]
  
  sobol_result <- sobolmartinez(model = R0_main_sobol, 
                                X1 = X1, 
                                X2 = X2, 
                                nboot = 100)
  
  sobols_first_Cres[i,-1] <- sobol_result$S$original
  sobols_total_Cres[i,-1] <- sobol_result$T$original
}
sobols_first_Cres[sobols_first_Cres < 0] <- 0
sobols_total_Cres[sobols_total_Cres < 0] <- 0

sobols_first_long_Cres <- sobols_first_Cres %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cres_sobol_first <- ggplot(sobols_first_long_Cres) +
  geom_line(mapping = aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) + 
  labs(x = "Temperature (°C)",
       title = "First-order Sobol' indices",
       color = "Trait") +
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

sobols_total_long_Cres <- sobols_total_Cres %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Cres_sobol_total <-ggplot(sobols_total_long_Cres) +
  geom_line(aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) +  
  labs(x = "Temperature (°C)",
       title = "Total-effect Sobol' indices",
       color = "Trait") + 
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Cres_sobol_first, plot_Cres_sobol_total, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. restuans", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/Sobol_Cres.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Cx. tarsalis

# calculate R0 samples from the trait samples
R0_main_Ctar = R0_main_f(biting_fit_Ctar, infprob_fit_pop, EIP_fit_Ctar, 
                         lf_fit_Ctar, omega, sur_fit_Ctar, ER, 
                         egg_viability_fit_pop, dev_fit_Ctar)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar <= 0) == ncol(R0_main_Ctar))){
  print("only zeros")
  R0_main_Ctar <- R0_main_Ctar[-which(rowSums(R0_main_Ctar <= 0) == ncol(R0_main_Ctar)),]  
}

# calculate statistics of relative R0 at each temperature
R0_main_Ctar_mean = apply(R0_main_Ctar, 2, mean)
R0_main_Ctar_0025 = apply(R0_main_Ctar, 2, quantile, probs=c(0.025))
R0_main_Ctar_0975 = apply(R0_main_Ctar, 2, quantile, probs=c(0.975))

# calculate Topt for each sample and calculate statistics
R0_main_Ctar_peaks = sapply(apply(R0_main_Ctar,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_peaks_mean = mean(R0_main_Ctar_peaks)
R0_main_Ctar_peaks_median = median(R0_main_Ctar_peaks)
R0_main_Ctar_peaks_0025 = quantile(R0_main_Ctar_peaks, probs=c(0.025))
R0_main_Ctar_peaks_0975 = quantile(R0_main_Ctar_peaks, probs=c(0.975))

# calculate Tmin for each sample and calculate statistics
R0_main_Ctar_Tmin = apply(R0_main_Ctar,1,FUN=Tmin)
R0_main_Ctar_Tmin_mean = mean(R0_main_Ctar_Tmin)
R0_main_Ctar_Tmin_median = median(R0_main_Ctar_Tmin)
R0_main_Ctar_Tmin_0025 = quantile(R0_main_Ctar_Tmin, probs=c(0.025))
R0_main_Ctar_Tmin_0975 = quantile(R0_main_Ctar_Tmin, probs=c(0.975))

# calculate Tmax for each sample and calculate statistics
R0_main_Ctar_Tmax = apply(R0_main_Ctar,1,FUN=Tmax)
R0_main_Ctar_Tmax_mean = mean(R0_main_Ctar_Tmax)
R0_main_Ctar_Tmax_median = median(R0_main_Ctar_Tmax)
R0_main_Ctar_Tmax_0025 = quantile(R0_main_Ctar_Tmax, probs=c(0.025))
R0_main_Ctar_Tmax_0975 = quantile(R0_main_Ctar_Tmax, probs=c(0.975))

# calculate R0 samples again but leaving one trait constant

# Biting rate

R0_main_Ctar_biting_const = R0_main_f(biting_data_mean, infprob_fit_pop, EIP_fit_Ctar, 
                                      lf_fit_Ctar, omega, sur_fit_Ctar, ER, 
                                      egg_viability_fit_pop, dev_fit_Ctar)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar_biting_const <= 0) == ncol(R0_main_Ctar_biting_const))){
  print("only zeros")
  R0_main_Ctar_biting_const <- R0_main_Ctar_biting_const[-which(rowSums(R0_main_Ctar_biting_const <= 0) == ncol(R0_main_Ctar_biting_const)),]  
}

R0_main_Ctar_biting_const_mean = apply(R0_main_Ctar_biting_const, 2, mean)

R0_main_Ctar_biting_const_peaks = sapply(apply(R0_main_Ctar_biting_const,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_biting_const_peaks_mean = mean(R0_main_Ctar_biting_const_peaks)
R0_main_Ctar_biting_const_peaks_0025 = quantile(R0_main_Ctar_biting_const_peaks, probs=c(0.025))
R0_main_Ctar_biting_const_peaks_0975 = quantile(R0_main_Ctar_biting_const_peaks, probs=c(0.975))

R0_main_Ctar_biting_const_Tmax = apply(R0_main_Ctar_biting_const,1,FUN=Tmax)
R0_main_Ctar_biting_const_Tmax_mean = mean(R0_main_Ctar_biting_const_Tmax)
R0_main_Ctar_biting_const_Tmax_0025 = quantile(R0_main_Ctar_biting_const_Tmax, probs=c(0.025))
R0_main_Ctar_biting_const_Tmax_0975 = quantile(R0_main_Ctar_biting_const_Tmax, probs=c(0.975))

R0_main_Ctar_biting_const_Tmin = apply(R0_main_Ctar_biting_const,1,FUN=Tmin)
R0_main_Ctar_biting_const_Tmin_mean = mean(R0_main_Ctar_biting_const_Tmin)
R0_main_Ctar_biting_const_Tmin_0025 = quantile(R0_main_Ctar_biting_const_Tmin, probs=c(0.025))
R0_main_Ctar_biting_const_Tmin_0975 = quantile(R0_main_Ctar_biting_const_Tmin, probs=c(0.975))

# Infection probability

R0_main_Ctar_infprob_const = R0_main_f(biting_fit_Ctar, infprob_data_mean, EIP_fit_Ctar, 
                                       lf_fit_Ctar, omega, sur_fit_Ctar, ER, 
                                       egg_viability_fit_pop, dev_fit_Ctar)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar_infprob_const <= 0) == ncol(R0_main_Ctar_infprob_const))){
  print("only zeros")
  R0_main_Ctar_infprob_const <- R0_main_Ctar_infprob_const[-which(rowSums(R0_main_Ctar_infprob_const <= 0) == ncol(R0_main_Ctar_infprob_const)),]  
}

R0_main_Ctar_infprob_const_mean = apply(R0_main_Ctar_infprob_const, 2, mean)

R0_main_Ctar_infprob_const_peaks = sapply(apply(R0_main_Ctar_infprob_const,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_infprob_const_peaks_mean = mean(R0_main_Ctar_infprob_const_peaks)
R0_main_Ctar_infprob_const_peaks_0025 = quantile(R0_main_Ctar_infprob_const_peaks, probs=c(0.025))
R0_main_Ctar_infprob_const_peaks_0975 = quantile(R0_main_Ctar_infprob_const_peaks, probs=c(0.975))

R0_main_Ctar_infprob_const_Tmax = apply(R0_main_Ctar_infprob_const,1,FUN=Tmax)
R0_main_Ctar_infprob_const_Tmax_mean = mean(R0_main_Ctar_infprob_const_Tmax)
R0_main_Ctar_infprob_const_Tmax_0025 = quantile(R0_main_Ctar_infprob_const_Tmax, probs=c(0.025))
R0_main_Ctar_infprob_const_Tmax_0975 = quantile(R0_main_Ctar_infprob_const_Tmax, probs=c(0.975))

R0_main_Ctar_infprob_const_Tmin = apply(R0_main_Ctar_infprob_const,1,FUN=Tmin)
R0_main_Ctar_infprob_const_Tmin_mean = mean(R0_main_Ctar_infprob_const_Tmin)
R0_main_Ctar_infprob_const_Tmin_0025 = quantile(R0_main_Ctar_infprob_const_Tmin, probs=c(0.025))
R0_main_Ctar_infprob_const_Tmin_0975 = quantile(R0_main_Ctar_infprob_const_Tmin, probs=c(0.975))

# EIP

R0_main_Ctar_EIP_const = R0_main_f(biting_fit_Ctar, infprob_fit_pop, EIP_mean, 
                                   lf_fit_Ctar, omega, sur_fit_Ctar, ER, 
                                   egg_viability_fit_pop, dev_fit_Ctar)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar_EIP_const <= 0) == ncol(R0_main_Ctar_EIP_const))){
  print("only zeros")
  R0_main_Ctar_EIP_const <- R0_main_Ctar_EIP_const[-which(rowSums(R0_main_Ctar_EIP_const <= 0) == ncol(R0_main_Ctar_EIP_const)),]  
}

R0_main_Ctar_EIP_const_mean = apply(R0_main_Ctar_EIP_const, 2, mean)

R0_main_Ctar_EIP_const_peaks = sapply(apply(R0_main_Ctar_EIP_const,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_EIP_const_peaks_mean = mean(R0_main_Ctar_EIP_const_peaks)
R0_main_Ctar_EIP_const_peaks_0025 = quantile(R0_main_Ctar_EIP_const_peaks, probs=c(0.025))
R0_main_Ctar_EIP_const_peaks_0975 = quantile(R0_main_Ctar_EIP_const_peaks, probs=c(0.975))

R0_main_Ctar_EIP_const_Tmax = apply(R0_main_Ctar_EIP_const,1,FUN=Tmax)
R0_main_Ctar_EIP_const_Tmax_mean = mean(R0_main_Ctar_EIP_const_Tmax)
R0_main_Ctar_EIP_const_Tmax_0025 = quantile(R0_main_Ctar_EIP_const_Tmax, probs=c(0.025))
R0_main_Ctar_EIP_const_Tmax_0975 = quantile(R0_main_Ctar_EIP_const_Tmax, probs=c(0.975))

R0_main_Ctar_EIP_const_Tmin = apply(R0_main_Ctar_EIP_const,1,FUN=Tmin)
R0_main_Ctar_EIP_const_Tmin_mean = mean(R0_main_Ctar_EIP_const_Tmin)
R0_main_Ctar_EIP_const_Tmin_0025 = quantile(R0_main_Ctar_EIP_const_Tmin, probs=c(0.025))
R0_main_Ctar_EIP_const_Tmin_0975 = quantile(R0_main_Ctar_EIP_const_Tmin, probs=c(0.975))

# Adult lifespan

R0_main_Ctar_lf_const = R0_main_f(biting_fit_Ctar, infprob_fit_pop, EIP_fit_Ctar, 
                                  lf_data_mean, omega, sur_fit_Ctar, ER, 
                                  egg_viability_fit_pop, dev_fit_Ctar)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar_lf_const <= 0) == ncol(R0_main_Ctar_lf_const))){
  print("only zeros")
  R0_main_Ctar_lf_const <- R0_main_Ctar_lf_const[-which(rowSums(R0_main_Ctar_lf_const <= 0) == ncol(R0_main_Ctar_lf_const)),]  
}

R0_main_Ctar_lf_const_mean = apply(R0_main_Ctar_lf_const, 2, mean)

R0_main_Ctar_lf_const_peaks = sapply(apply(R0_main_Ctar_lf_const,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_lf_const_peaks_mean = mean(R0_main_Ctar_lf_const_peaks)
R0_main_Ctar_lf_const_peaks_0025 = quantile(R0_main_Ctar_lf_const_peaks, probs=c(0.025))
R0_main_Ctar_lf_const_peaks_0975 = quantile(R0_main_Ctar_lf_const_peaks, probs=c(0.975))

R0_main_Ctar_lf_const_Tmax = apply(R0_main_Ctar_lf_const,1,FUN=Tmax)
R0_main_Ctar_lf_const_Tmax_mean = mean(R0_main_Ctar_lf_const_Tmax)
R0_main_Ctar_lf_const_Tmax_0025 = quantile(R0_main_Ctar_lf_const_Tmax, probs=c(0.025))
R0_main_Ctar_lf_const_Tmax_0975 = quantile(R0_main_Ctar_lf_const_Tmax, probs=c(0.975))

R0_main_Ctar_lf_const_Tmin = apply(R0_main_Ctar_lf_const,1,FUN=Tmin)
R0_main_Ctar_lf_const_Tmin_mean = mean(R0_main_Ctar_lf_const_Tmin)
R0_main_Ctar_lf_const_Tmin_0025 = quantile(R0_main_Ctar_lf_const_Tmin, probs=c(0.025))
R0_main_Ctar_lf_const_Tmin_0975 = quantile(R0_main_Ctar_lf_const_Tmin, probs=c(0.975))

# Juvenile survival

R0_main_Ctar_sur_const = R0_main_f(biting_fit_Ctar, infprob_fit_pop, EIP_fit_Ctar, 
                                   lf_fit_Ctar, omega, sur_data_mean, ER, 
                                   egg_viability_fit_pop, dev_fit_Ctar)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar_sur_const <= 0) == ncol(R0_main_Ctar_sur_const))){
  print("only zeros")
  R0_main_Ctar_sur_const <- R0_main_Ctar_sur_const[-which(rowSums(R0_main_Ctar_sur_const <= 0) == ncol(R0_main_Ctar_sur_const)),]  
}

R0_main_Ctar_sur_const_mean = apply(R0_main_Ctar_sur_const, 2, mean)

R0_main_Ctar_sur_const_peaks = sapply(apply(R0_main_Ctar_sur_const,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_sur_const_peaks_mean = mean(R0_main_Ctar_sur_const_peaks)
R0_main_Ctar_sur_const_peaks_0025 = quantile(R0_main_Ctar_sur_const_peaks, probs=c(0.025))
R0_main_Ctar_sur_const_peaks_0975 = quantile(R0_main_Ctar_sur_const_peaks, probs=c(0.975))

R0_main_Ctar_sur_const_Tmax = apply(R0_main_Ctar_sur_const,1,FUN=Tmax)
R0_main_Ctar_sur_const_Tmax_mean = mean(R0_main_Ctar_sur_const_Tmax)
R0_main_Ctar_sur_const_Tmax_0025 = quantile(R0_main_Ctar_sur_const_Tmax, probs=c(0.025))
R0_main_Ctar_sur_const_Tmax_0975 = quantile(R0_main_Ctar_sur_const_Tmax, probs=c(0.975))

R0_main_Ctar_sur_const_Tmin = apply(R0_main_Ctar_sur_const,1,FUN=Tmin)
R0_main_Ctar_sur_const_Tmin_mean = mean(R0_main_Ctar_sur_const_Tmin)
R0_main_Ctar_sur_const_Tmin_0025 = quantile(R0_main_Ctar_sur_const_Tmin, probs=c(0.025))
R0_main_Ctar_sur_const_Tmin_0975 = quantile(R0_main_Ctar_sur_const_Tmin, probs=c(0.975))

# Egg viability

R0_main_Ctar_EV_const = R0_main_f(biting_fit_Ctar, infprob_fit_pop, EIP_fit_Ctar, 
                                  lf_fit_Ctar, omega, sur_fit_Ctar, ER, 
                                  egg_viability_data_mean, dev_fit_Ctar)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar_EV_const <= 0) == ncol(R0_main_Ctar_EV_const))){
  print("only zeros")
  R0_main_Ctar_EV_const <- R0_main_Ctar_EV_const[-which(rowSums(R0_main_Ctar_EV_const <= 0) == ncol(R0_main_Ctar_EV_const)),]  
}

R0_main_Ctar_EV_const_mean = apply(R0_main_Ctar_EV_const, 2, mean)

R0_main_Ctar_EV_const_peaks = sapply(apply(R0_main_Ctar_EV_const,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_EV_const_peaks_mean = mean(R0_main_Ctar_EV_const_peaks)
R0_main_Ctar_EV_const_peaks_0025 = quantile(R0_main_Ctar_EV_const_peaks, probs=c(0.025))
R0_main_Ctar_EV_const_peaks_0975 = quantile(R0_main_Ctar_EV_const_peaks, probs=c(0.975))

R0_main_Ctar_EV_const_Tmax = apply(R0_main_Ctar_EV_const,1,FUN=Tmax)
R0_main_Ctar_EV_const_Tmax_mean = mean(R0_main_Ctar_EV_const_Tmax)
R0_main_Ctar_EV_const_Tmax_0025 = quantile(R0_main_Ctar_EV_const_Tmax, probs=c(0.025))
R0_main_Ctar_EV_const_Tmax_0975 = quantile(R0_main_Ctar_EV_const_Tmax, probs=c(0.975))

R0_main_Ctar_EV_const_Tmin = apply(R0_main_Ctar_EV_const,1,FUN=Tmin)
R0_main_Ctar_EV_const_Tmin_mean = mean(R0_main_Ctar_EV_const_Tmin)
R0_main_Ctar_EV_const_Tmin_0025 = quantile(R0_main_Ctar_EV_const_Tmin, probs=c(0.025))
R0_main_Ctar_EV_const_Tmin_0975 = quantile(R0_main_Ctar_EV_const_Tmin, probs=c(0.975))

# Juvenile development

R0_main_Ctar_dev_const = R0_main_f(biting_fit_Ctar, infprob_fit_pop, EIP_fit_Ctar, 
                                   lf_fit_Ctar, omega, sur_fit_Ctar, ER, 
                                   egg_viability_fit_pop, dev_data_mean)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar_dev_const <= 0) == ncol(R0_main_Ctar_dev_const))){
  print("only zeros")
  R0_main_Ctar_dev_const <- R0_main_Ctar_dev_const[-which(rowSums(R0_main_Ctar_dev_const <= 0) == ncol(R0_main_Ctar_dev_const)),]  
}

R0_main_Ctar_dev_const_mean = apply(R0_main_Ctar_dev_const, 2, mean)

R0_main_Ctar_dev_const_peaks = sapply(apply(R0_main_Ctar_dev_const,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_dev_const_peaks_mean = mean(R0_main_Ctar_dev_const_peaks)
R0_main_Ctar_dev_const_peaks_0025 = quantile(R0_main_Ctar_dev_const_peaks, probs=c(0.025))
R0_main_Ctar_dev_const_peaks_0975 = quantile(R0_main_Ctar_dev_const_peaks, probs=c(0.975))

R0_main_Ctar_dev_const_Tmax = apply(R0_main_Ctar_dev_const,1,FUN=Tmax)
R0_main_Ctar_dev_const_Tmax_mean = mean(R0_main_Ctar_dev_const_Tmax)
R0_main_Ctar_dev_const_Tmax_0025 = quantile(R0_main_Ctar_dev_const_Tmax, probs=c(0.025))
R0_main_Ctar_dev_const_Tmax_0975 = quantile(R0_main_Ctar_dev_const_Tmax, probs=c(0.975))

R0_main_Ctar_dev_const_Tmin = apply(R0_main_Ctar_dev_const,1,FUN=Tmin)
R0_main_Ctar_dev_const_Tmin_mean = mean(R0_main_Ctar_dev_const_Tmin)
R0_main_Ctar_dev_const_Tmin_0025 = quantile(R0_main_Ctar_dev_const_Tmin, probs=c(0.025))
R0_main_Ctar_dev_const_Tmin_0975 = quantile(R0_main_Ctar_dev_const_Tmin, probs=c(0.975))

# Summarize results for Cx. pip. molestus in dataframe

df_R0_main_Ctar_const <- data.frame(x = temp,
                                    y = c(R0_main_Ctar_biting_const_mean/max(R0_main_Ctar_biting_const_mean),
                                          R0_main_Ctar_infprob_const_mean/max(R0_main_Ctar_infprob_const_mean),
                                          R0_main_Ctar_EIP_const_mean/max(R0_main_Ctar_EIP_const_mean),
                                          R0_main_Ctar_lf_const_mean/max(R0_main_Ctar_lf_const_mean),
                                          R0_main_Ctar_sur_const_mean/max(R0_main_Ctar_sur_const_mean),
                                          R0_main_Ctar_EV_const_mean/max(R0_main_Ctar_EV_const_mean),
                                          R0_main_Ctar_dev_const_mean/max(R0_main_Ctar_dev_const_mean),
                                          R0_main_Ctar_mean/max(R0_main_Ctar_mean)
                                    ),
                                    trait = c(rep("3Biting rate", length(temp)),
                                              rep("2Inf. prob.", length(temp)),
                                              rep("1EIP", length(temp)),
                                              rep("5Lifespan", length(temp)),
                                              rep("6Juv. survival", length(temp)),
                                              rep("4Egg viab.", length(temp)),
                                              rep("7Juv. dev. rate", length(temp)),
                                              rep("8None", length(temp))
                                    )
)

df_R0_main_Ctar_const_stats <- data.frame(peaks = c(R0_main_Ctar_biting_const_peaks_mean,
                                                    R0_main_Ctar_infprob_const_peaks_mean,
                                                    R0_main_Ctar_EIP_const_peaks_mean,
                                                    R0_main_Ctar_lf_const_peaks_mean,
                                                    R0_main_Ctar_sur_const_peaks_mean,
                                                    R0_main_Ctar_EV_const_peaks_mean,
                                                    R0_main_Ctar_dev_const_peaks_mean,
                                                    R0_main_Ctar_peaks_mean),
                                          peaks_lower = c(R0_main_Ctar_biting_const_peaks_0025,
                                                          R0_main_Ctar_infprob_const_peaks_0025,
                                                          R0_main_Ctar_EIP_const_peaks_0025,
                                                          R0_main_Ctar_lf_const_peaks_0025,
                                                          R0_main_Ctar_sur_const_peaks_0025,
                                                          R0_main_Ctar_EV_const_peaks_0025,
                                                          R0_main_Ctar_dev_const_peaks_0025,
                                                          R0_main_Ctar_peaks_0025),
                                          peaks_upper = c(R0_main_Ctar_biting_const_peaks_0975,
                                                          R0_main_Ctar_infprob_const_peaks_0975,
                                                          R0_main_Ctar_EIP_const_peaks_0975,
                                                          R0_main_Ctar_lf_const_peaks_0975,
                                                          R0_main_Ctar_sur_const_peaks_0975,
                                                          R0_main_Ctar_EV_const_peaks_0975,
                                                          R0_main_Ctar_dev_const_peaks_0975,
                                                          R0_main_Ctar_peaks_0975),
                                          Tmin = c(R0_main_Ctar_biting_const_Tmin_mean,
                                                   R0_main_Ctar_infprob_const_Tmin_mean,
                                                   R0_main_Ctar_EIP_const_Tmin_mean,
                                                   R0_main_Ctar_lf_const_Tmin_mean,
                                                   R0_main_Ctar_sur_const_Tmin_mean,
                                                   R0_main_Ctar_EV_const_Tmin_mean,
                                                   R0_main_Ctar_dev_const_Tmin_mean,
                                                   R0_main_Ctar_Tmin_mean),
                                          Tmin_lower = c(R0_main_Ctar_biting_const_Tmin_0025,
                                                         R0_main_Ctar_infprob_const_Tmin_0025,
                                                         R0_main_Ctar_EIP_const_Tmin_0025,
                                                         R0_main_Ctar_lf_const_Tmin_0025,
                                                         R0_main_Ctar_sur_const_Tmin_0025,
                                                         R0_main_Ctar_EV_const_Tmin_0025,
                                                         R0_main_Ctar_dev_const_Tmin_0025,
                                                         R0_main_Ctar_Tmin_0025),
                                          Tmin_upper = c(R0_main_Ctar_biting_const_Tmin_0975,
                                                         R0_main_Ctar_infprob_const_Tmin_0975,
                                                         R0_main_Ctar_EIP_const_Tmin_0975,
                                                         R0_main_Ctar_lf_const_Tmin_0975,
                                                         R0_main_Ctar_sur_const_Tmin_0975,
                                                         R0_main_Ctar_EV_const_Tmin_0975,
                                                         R0_main_Ctar_dev_const_Tmin_0975,
                                                         R0_main_Ctar_Tmin_0975),
                                          Tmax = c(R0_main_Ctar_biting_const_Tmax_mean,
                                                   R0_main_Ctar_infprob_const_Tmax_mean,
                                                   R0_main_Ctar_EIP_const_Tmax_mean,
                                                   R0_main_Ctar_lf_const_Tmax_mean,
                                                   R0_main_Ctar_sur_const_Tmax_mean,
                                                   R0_main_Ctar_EV_const_Tmax_mean,
                                                   R0_main_Ctar_dev_const_Tmax_mean,
                                                   R0_main_Ctar_Tmax_mean),
                                          Tmax_lower = c(R0_main_Ctar_biting_const_Tmax_0025,
                                                         R0_main_Ctar_infprob_const_Tmax_0025,
                                                         R0_main_Ctar_EIP_const_Tmax_0025,
                                                         R0_main_Ctar_lf_const_Tmax_0025,
                                                         R0_main_Ctar_sur_const_Tmax_0025,
                                                         R0_main_Ctar_EV_const_Tmax_0025,
                                                         R0_main_Ctar_dev_const_Tmax_0025,
                                                         R0_main_Ctar_Tmax_0025),
                                          Tmax_upper = c(R0_main_Ctar_biting_const_Tmax_0975,
                                                         R0_main_Ctar_infprob_const_Tmax_0975,
                                                         R0_main_Ctar_EIP_const_Tmax_0975,
                                                         R0_main_Ctar_lf_const_Tmax_0975,
                                                         R0_main_Ctar_sur_const_Tmax_0975,
                                                         R0_main_Ctar_EV_const_Tmax_0975,
                                                         R0_main_Ctar_dev_const_Tmax_0975,
                                                         R0_main_Ctar_Tmax_0975),
                                          trait = c("3Biting rate", "2Inf. prob.", "1EIP",
                                                    "5Lifespan", "6Juv. survival", "4Egg viab.", "7Juv. dev. rate",
                                                    "8None")
)

plot_Ctar_const <- ggplot() +
  geom_line(df_R0_main_Ctar_const, mapping = aes(x = x, y = y, color = trait), linewidth=0.6) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. restuans")))) +
  labs(x = "Temperature (°C)",
       color = "Constant trait", 
       title = expression("Mean temperature response of " * R[0]^rel)) +
  scale_color_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                                  "Lifespan", "Juv. survival", "Juv. dev. rate",
                                  "None")) +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_Ctar_stats <- ggplot(df_R0_main_Ctar_const_stats) +
  geom_pointrange(aes(x=peaks, xmin=peaks_lower, xmax=peaks_upper, y=trait, color=trait), linewidth=1.2, alpha=0.5, orientation = "y") +
  geom_pointrange(aes(x=Tmin, xmin=Tmin_lower, xmax=Tmin_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  geom_pointrange(aes(x=Tmax, xmin=Tmax_lower, xmax=Tmax_upper, y = trait, color=trait), linewidth=0.6, orientation = "y") +
  scale_x_continuous(breaks = seq(0,45,5), limits = c(2,41)) +
  guides(color = "none") +
  theme_bw() +
  labs(x = "Temperature (°C)",
       title = "      Temperature limits and optimal temperature") +
  scale_y_discrete(labels = c("EIP", "Mosq. inf. prob.", "Biting rate", "Egg viab.",
                              "Lifespan", "Juv. survival", "Juv. dev. rate",
                              "None")) +
  theme(plot.margin = unit(c(0.4, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Ctar_const, plot_Ctar_stats, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. tarsalis", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/R0_const_Ctar.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")

# Sobol' analysis

r = 151:301
n = length(r)

sobols_first_Ctar <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

sobols_total_Ctar <- data.frame(x = temp[r],
                                a = rep(NA, n),
                                bm = rep(NA, n),
                                EIP = rep(NA, n),
                                lf = rep(NA, n),
                                surJ = rep(NA, n),
                                EV = rep(NA, n),
                                devJ = rep(NA, n))

for(i in 1:n){
  input_matrix <- data.frame(a = biting_fit_Ctar[,r[i]], bM = infprob_fit_pop[,r[i]], EIP = EIP_fit_Ctar[,r[i]], 
                             lf = lf_fit_Ctar[,r[i]], surJ = sur_fit_Ctar[,r[i]], 
                             EV = egg_viability_fit_pop[,r[i]], devJ = dev_fit_Ctar[,r[i]])
  
  X1 <- input_matrix[X_index1, ]
  X2 <- input_matrix[X_index2, ]
  
  sobol_result <- sobolmartinez(model = R0_main_sobol, 
                                X1 = X1, 
                                X2 = X2, 
                                nboot = 100)
  
  sobols_first_Ctar[i,-1] <- sobol_result$S$original
  sobols_total_Ctar[i,-1] <- sobol_result$T$original
}
sobols_first_Ctar[sobols_first_Ctar < 0] <- 0
sobols_total_Ctar[sobols_total_Ctar < 0] <- 0

sobols_first_long_Ctar <- sobols_first_Ctar %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Ctar_sobol_first <- ggplot(sobols_first_long_Ctar) +
  geom_line(mapping = aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) + 
  labs(x = "Temperature (°C)",
       title = "First-order Sobol' indices",
       color = "Trait") +
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

sobols_total_long_Ctar <- sobols_total_Ctar %>%
  pivot_longer(cols = -x,               
               names_to = "param",   
               values_to = "y") 

plot_Ctar_sobol_total <-ggplot(sobols_total_long_Ctar) +
  geom_line(aes(x = x, y = y, color = param)) +
  ylim(0, 1) + 
  xlim(15,30) +  
  labs(x = "Temperature (°C)",
       title = "Total-effect Sobol' indices",
       color = "Trait") + 
  scale_color_discrete(labels = c("Biting rate", "Mosq. inf. prob.", "Juv. dev. rate", "EIP",
                                  "Egg viab.", "Lifespan", "Juv. survival")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),  
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.text = element_text(size=10),
        plot.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot <- ggarrange(plot_Ctar_sobol_first, plot_Ctar_sobol_total, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="bottom", labels=c("A", "B"), font.label = list(size = 12))

plot <- annotate_figure(plot,
                        top = text_grob("Cx. tarsalis", 
                                        size = 11, face = "bold.italic"))

plot

#ggsave("Figures/Sobol_Ctar.tiff", 
#       plot = plot, 
#       width = 6.3, height = 4, units = "in", dpi = 600,
#       compression = "lzw")
