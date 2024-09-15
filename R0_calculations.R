### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to calculate and visualize temperature responses of 
### the relative basic reproduction number of West Nile virus in different
### Culex mosquitoes. The calculations are based on temperature responses
### of different mosquito-pathogen traits that were derived with a Bayesian
### hierarchical model implemented in STAN. 
###
### The script includes the necessary code to reproduce Figures 3, 4, SI4.5
###

# load libraries
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

## load fitted models, extract trait temperature response samples for each species
## and population-level estimates (which we use to substitute the temperature 
## response of a species when it lacks data for a specific trait), and remove models
## afterwards to save local memory

# EIP
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
biting_fit <- readRDS("model_fits/biting_rate_fit.rds")
biting_fit_Cpal <- rstan::extract(biting_fit, permuted=T)$f_new[,2,]
biting_fit_Cpip <- rstan::extract(biting_fit, permuted=T)$f_new[,3,]
biting_fit_Cqui <- rstan::extract(biting_fit, permuted=T)$f_new[,4,]
biting_fit_Ctar <- rstan::extract(biting_fit, permuted=T)$f_new[,5,]
biting_fit_pop <- rstan::extract(biting_fit, permuted=T)$f_new_spec
rm(biting_fit)

# Adult mosquito lifespan
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
dev_egg_fit <- readRDS("model_fits/egg_development_fit.rds")
dev_egg_fit_Cmol <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,2,]
dev_egg_fit_Cpal <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,3,]
dev_egg_fit_Cpip <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,4,]
dev_egg_fit_Cqui <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,5,]
dev_egg_fit_Cres <- rstan::extract(dev_egg_fit, permuted=T)$f_new[,6,]
rm(dev_egg_fit)

# Juvenile development rate
dev_fit <- readRDS("model_fits/juvenile_development_fit.rds")
dev_fit_Cmol <- rstan::extract(dev_fit, permuted=T)$f_new[,5,]
dev_fit_Cpal <- rstan::extract(dev_fit, permuted=T)$f_new[,6,]
dev_fit_Cpip <- rstan::extract(dev_fit, permuted=T)$f_new[,7,]
dev_fit_Cqui <- rstan::extract(dev_fit, permuted=T)$f_new[,8,]
dev_fit_Cres <- rstan::extract(dev_fit, permuted=T)$f_new[,9,]
dev_fit_Ctar <- rstan::extract(dev_fit, permuted=T)$f_new[,11,]
rm(dev_fit)

# Egg viability
egg_viability_fit <- readRDS("model_fits/egg_viability_fit.rds")
egg_viability_fit_Cmol <- rstan::extract(egg_viability_fit, permuted=T)$f_new[,1,]
egg_viability_fit_Cpal <- rstan::extract(egg_viability_fit, permuted=T)$f_new[,2,]
egg_viability_fit_Cqui <- rstan::extract(egg_viability_fit, permuted=T)$f_new[,3,]
egg_viability_fit_pop <- rstan::extract(egg_viability_fit, permuted=T)$f_new_spec
rm(egg_viability_fit)

# Mosquito infection probability
infprob_fit <- readRDS("model_fits/infection_probability_fit.rds")
infprob_fit_pop <- rstan::extract(infprob_fit, permuted=T)$f_new_spec
rm(infprob_fit)
infprob_fit_reduced <- readRDS("model_fits/infection_probability_fit_reduced_experiment_variability.rds")
infprob_fit_pop_red <- rstan::extract(infprob_fit_reduced, permuted=T)$f_new_spec
rm(infprob_fit_reduced)
infprob_fit_increased <- readRDS("model_fits/infection_probability_fit_increased_experiment_variability.rds")
infprob_fit_pop_inc <- rstan::extract(infprob_fit_increased, permuted=T)$f_new_spec
rm(infprob_fit_increased)

# Juvenile survival
sur_fit <- readRDS("model_fits/juvenile_survival_fit.rds")
sur_fit_Cmol <- rstan::extract(sur_fit, permuted=T)$f_new[,7,]
sur_fit_Cpal <- rstan::extract(sur_fit, permuted=T)$f_new[,8,]
sur_fit_Cpip <- rstan::extract(sur_fit, permuted=T)$f_new[,9,]
sur_fit_Cqui <- rstan::extract(sur_fit, permuted=T)$f_new[,10,]
sur_fit_Cres <- rstan::extract(sur_fit, permuted=T)$f_new[,11,]
sur_fit_Ctar <- rstan::extract(sur_fit, permuted=T)$f_new[,13,]
rm(sur_fit)

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

## relative R0 calculations for each species including calculation of statistics 
## for Tmin, Tmax, and Topt (temperature where relative R0 peaks) + plotting

# Cx. pipiens

# calculate R0 samples from the trait samples
R0_main_Cpip = R0_main_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, lf_fit_Cpip, omega, sur_fit_Cpip, ER, egg_viability_fit_pop, dev_fit_Cpip)

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

# collect information in dataframe and normalize by value at optimal temperature
df_Cpip <- data.frame(x = temp, mean = c(R0_main_Cpip_mean/max(R0_main_Cpip_mean)),
                      lowerCI = c(R0_main_Cpip_0025/max(R0_main_Cpip_mean)),
                      upperCI = c(R0_main_Cpip_0975/max(R0_main_Cpip_mean)),
                      Species = c(rep("pip.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_main_Cpip_peaks = sapply(apply(R0_main_Cpip,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_peaks_mean = mean(R0_main_Cpip_peaks)
R0_main_Cpip_peaks_median = median(R0_main_Cpip_peaks)
R0_main_Cpip_peaks_0025 = quantile(R0_main_Cpip_peaks, probs=c(0.025))
R0_main_Cpip_peaks_0975 = quantile(R0_main_Cpip_peaks, probs=c(0.975))

# check if any sample has a positive value at lowest temperature which would make 
# Tmin calculations nonsensical
any(R0_main_Cpip[,1]>0)

# calculate Tmin for each sample and calculate statistics
R0_main_Cpip_Tmin = apply(R0_main_Cpip,1,FUN=Tmin)
R0_main_Cpip_Tmin_mean = mean(R0_main_Cpip_Tmin)
R0_main_Cpip_Tmin_median =median(R0_main_Cpip_Tmin)
R0_main_Cpip_Tmin_0025 = quantile(R0_main_Cpip_Tmin, probs=c(0.025))
R0_main_Cpip_Tmin_0975 = quantile(R0_main_Cpip_Tmin, probs=c(0.975))

# check if any sample has a positive value at highest temperature which would make 
# Tmax calculations nonsensical
any(R0_main_Cpip[,451]>0)

# calculate Tmax for each sample and calculate statistics
R0_main_Cpip_Tmax = apply(R0_main_Cpip,1,FUN=Tmax)
R0_main_Cpip_Tmax_mean = mean(R0_main_Cpip_Tmax)
R0_main_Cpip_Tmax_median =median(R0_main_Cpip_Tmax)
R0_main_Cpip_Tmax_0025 = quantile(R0_main_Cpip_Tmax, probs=c(0.025))
R0_main_Cpip_Tmax_0975 = quantile(R0_main_Cpip_Tmax, probs=c(0.975))

# plot mean relative R0 temperature response with mean and 95% CIs for 
# Tmin, Topt, Tmax
plot_Cpip <- ggplot() +
  geom_line(df_Cpip, mapping = aes(x = x, y = mean),linewidth=0.8, color="red") + 
  geom_line(data = data.frame(x = c(R0_main_Cpip_peaks_0025, R0_main_Cpip_peaks_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cpip_Tmin_0025, R0_main_Cpip_Tmin_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cpip_Tmax_0025, R0_main_Cpip_Tmax_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_point(data = data.frame(x = c(R0_main_Cpip_Tmin_mean, R0_main_Cpip_peaks_mean, R0_main_Cpip_Tmax_mean), y = 0), aes(x = x, y = y), size=2) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(5,40)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. pipiens")))) +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_blank(), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


# Cx. quinquefasciatus

# calculate R0 samples from the trait samples
R0_main_Cqui = R0_main_f(biting_fit_Cqui, infprob_fit_pop, EIP_fit_pop, lf_fit_Cqui, omega, sur_fit_Cqui, ER, egg_viability_fit_Cqui, dev_fit_Cqui)

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

# collect information in dataframe and normalize by value at optimal temperature
df_Cqui <- data.frame(x = temp, mean = c(R0_main_Cqui_mean/max(R0_main_Cqui_mean)),
                      lowerCI = c(R0_main_Cqui_0025/max(R0_main_Cqui_mean)),
                      upperCI = c(R0_main_Cqui_0975/max(R0_main_Cqui_mean)),
                      Species = c(rep("qui.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_main_Cqui_peaks = sapply(apply(R0_main_Cqui,1,which.max), FUN=function(x)temp[x])
R0_main_Cqui_peaks_mean = mean(R0_main_Cqui_peaks)
R0_main_Cqui_peaks_median =median(R0_main_Cqui_peaks)
R0_main_Cqui_peaks_0025 = quantile(R0_main_Cqui_peaks, probs=c(0.025))
R0_main_Cqui_peaks_0975 = quantile(R0_main_Cqui_peaks, probs=c(0.975))

# check if any sample has a positive value at lowest temperature which would make 
# Tmin calculations nonsensical
any(R0_main_Cqui[,1]>0)

# calculate Tmin for each sample and calculate statistics
R0_main_Cqui_Tmin = apply(R0_main_Cqui,1,FUN=Tmin)
R0_main_Cqui_Tmin_mean = mean(R0_main_Cqui_Tmin)
R0_main_Cqui_Tmin_median =median(R0_main_Cqui_Tmin)
R0_main_Cqui_Tmin_0025 = quantile(R0_main_Cqui_Tmin, probs=c(0.025))
R0_main_Cqui_Tmin_0975 = quantile(R0_main_Cqui_Tmin, probs=c(0.975))

# check if any sample has a positive value at highest temperature which would make 
# Tmax calculations nonsensical
any(R0_main_Cqui[,451]>0)

# calculate Tmax for each sample and calculate statistics
R0_main_Cqui_Tmax = apply(R0_main_Cqui,1,FUN=Tmax)
R0_main_Cqui_Tmax_mean = mean(R0_main_Cqui_Tmax)
R0_main_Cqui_Tmax_median =median(R0_main_Cqui_Tmax)
R0_main_Cqui_Tmax_0025 = quantile(R0_main_Cqui_Tmax, probs=c(0.025))
R0_main_Cqui_Tmax_0975 = quantile(R0_main_Cqui_Tmax, probs=c(0.975))

# plot mean relative R0 temperature response with mean and 95% CIs for 
# Tmin, Topt, Tmax
plot_Cqui <- ggplot() +
  geom_line(df_Cqui, mapping = aes(x = x, y = mean),linewidth=0.8, color="red") + 
  geom_line(data = data.frame(x = c(R0_main_Cqui_peaks_0025, R0_main_Cqui_peaks_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cqui_Tmin_0025, R0_main_Cqui_Tmin_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cqui_Tmax_0025, R0_main_Cqui_Tmax_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_point(data = data.frame(x = c(R0_main_Cqui_Tmin_mean, R0_main_Cqui_peaks_mean, R0_main_Cqui_Tmax_mean), y = 0), aes(x = x, y = y), size=2) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(5,40)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. quinquefasciatus")))) +
  theme(axis.text = element_text(size = 14),  
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

# Cx. pipiens molestus

# calculate R0 samples from the trait samples
R0_main_Cmol = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, lf_fit_Cmol, omega, sur_fit_Cmol, ER, egg_viability_fit_Cmol, dev_fit_Cmol)

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

# collect information in dataframe and normalize by value at optimal temperature
df_Cmol <- data.frame(x = temp, mean = c(R0_main_Cmol_mean/max(R0_main_Cmol_mean)),
                      lowerCI = c(R0_main_Cmol_0025/max(R0_main_Cmol_mean)),
                      upperCI = c(R0_main_Cmol_0975/max(R0_main_Cmol_mean)),
                      Species = c(rep("mol.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_main_Cmol_peaks = sapply(apply(R0_main_Cmol,1,which.max), FUN=function(x)temp[x])
R0_main_Cmol_peaks_mean = mean(R0_main_Cmol_peaks)
R0_main_Cmol_peaks_median =median(R0_main_Cmol_peaks)
R0_main_Cmol_peaks_0025 = quantile(R0_main_Cmol_peaks, probs=c(0.025))
R0_main_Cmol_peaks_0975 = quantile(R0_main_Cmol_peaks, probs=c(0.975))

# check if any sample has a positive value at lowest temperature which would make 
# Tmin calculations nonsensical
any(R0_main_Cmol[,1]>0)

# calculate Tmin for each sample and calculate statistics
R0_main_Cmol_Tmin = apply(R0_main_Cmol,1,FUN=Tmin)
R0_main_Cmol_Tmin_mean = mean(R0_main_Cmol_Tmin)
R0_main_Cmol_Tmin_median =median(R0_main_Cmol_Tmin)
R0_main_Cmol_Tmin_0025 = quantile(R0_main_Cmol_Tmin, probs=c(0.025))
R0_main_Cmol_Tmin_0975 = quantile(R0_main_Cmol_Tmin, probs=c(0.975))

# check if any sample has a positive value at highest temperature which would make 
# Tmax calculations nonsensical
any(R0_main_Cmol[,451]>0)

# calculate Tmax for each sample and calculate statistics
R0_main_Cmol_Tmax = apply(R0_main_Cmol,1,FUN=Tmax)
R0_main_Cmol_Tmax_mean = mean(R0_main_Cmol_Tmax)
R0_main_Cmol_Tmax_median =median(R0_main_Cmol_Tmax)
R0_main_Cmol_Tmax_0025 = quantile(R0_main_Cmol_Tmax, probs=c(0.025))
R0_main_Cmol_Tmax_0975 = quantile(R0_main_Cmol_Tmax, probs=c(0.975))

# plot mean relative R0 temperature response with mean and 95% CIs for 
# Tmin, Topt, Tmax
plot_Cmol <- ggplot() +
  geom_line(df_Cmol, mapping = aes(x = x, y = mean),linewidth=0.8, color="red") + 
  geom_line(data = data.frame(x = c(R0_main_Cmol_peaks_0025, R0_main_Cmol_peaks_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cmol_Tmin_0025, R0_main_Cmol_Tmin_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cmol_Tmax_0025, R0_main_Cmol_Tmax_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_point(data = data.frame(x = c(R0_main_Cmol_Tmin_mean, R0_main_Cmol_peaks_mean, R0_main_Cmol_Tmax_mean), y = 0), aes(x = x, y = y), size=2) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(5,40)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. pipiens molestus")))) +
  theme(axis.text = element_text(size = 14),  
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

# Cx. pipiens pallens

# calculate R0 samples from the trait samples
R0_main_Cpal = R0_main_f(biting_fit_Cpal, infprob_fit_pop, EIP_fit_pop, lf_fit_Cpal, omega, sur_fit_Cpal, ER, egg_viability_fit_Cpal, dev_fit_Cpal)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cpal <= 0) == ncol(R0_main_Cpal))){
  print("lol")
  R0_main_Cpal <- R0_main_Cpal[-which(rowSums(R0_main_Cpal <= 0) == ncol(R0_main_Cpal)),]
}

# calculate statistics of relative R0 at each temperature
R0_main_Cpal_mean = apply(R0_main_Cpal, 2, mean)
R0_main_Cpal_0025 = apply(R0_main_Cpal, 2, quantile, probs=c(0.025))
R0_main_Cpal_0975 = apply(R0_main_Cpal, 2, quantile, probs=c(0.975))

# collect information in dataframe and normalize by value at optimal temperature
df_Cpal <- data.frame(x = temp, mean = c(R0_main_Cpal_mean/max(R0_main_Cpal_mean)),
                      lowerCI = c(R0_main_Cpal_0025/max(R0_main_Cpal_mean)),
                      upperCI = c(R0_main_Cpal_0975/max(R0_main_Cpal_mean)),
                      Species = c(rep("pal.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_main_Cpal_peaks = sapply(apply(R0_main_Cpal,1,which.max), FUN=function(x)temp[x])
R0_main_Cpal_peaks_mean = mean(R0_main_Cpal_peaks)
R0_main_Cpal_peaks_median =median(R0_main_Cpal_peaks)
R0_main_Cpal_peaks_0025 = quantile(R0_main_Cpal_peaks, probs=c(0.025))
R0_main_Cpal_peaks_0975 = quantile(R0_main_Cpal_peaks, probs=c(0.975))

# check if any sample has a positive value at lowest temperature which would make 
# Tmin calculations nonsensical
any(R0_main_Cpal[,1]>0)

# calculate Tmin for each sample and calculate statistics
R0_main_Cpal_Tmin = apply(R0_main_Cpal,1,FUN=Tmin)
R0_main_Cpal_Tmin_mean = mean(R0_main_Cpal_Tmin)
R0_main_Cpal_Tmin_median = median(R0_main_Cpal_Tmin)
R0_main_Cpal_Tmin_0025 = quantile(R0_main_Cpal_Tmin, probs=c(0.025))
R0_main_Cpal_Tmin_0975 = quantile(R0_main_Cpal_Tmin, probs=c(0.975))

# check if any sample has a positive value at highest temperature which would make 
# Tmax calculations nonsensical
any(R0_main_Cpal[,451]>0)

# calculate Tmax for each sample and calculate statistics
R0_main_Cpal_Tmax = apply(R0_main_Cpal,1,FUN=Tmax)
R0_main_Cpal_Tmax_mean = mean(R0_main_Cpal_Tmax)
R0_main_Cpal_Tmax_median =median(R0_main_Cpal_Tmax)
R0_main_Cpal_Tmax_0025 = quantile(R0_main_Cpal_Tmax, probs=c(0.025))
R0_main_Cpal_Tmax_0975 = quantile(R0_main_Cpal_Tmax, probs=c(0.975))

# plot mean relative R0 temperature response with mean and 95% CIs for 
# Tmin, Topt, Tmax
plot_Cpal <- ggplot() +
  geom_line(df_Cpal, mapping = aes(x = x, y = mean),linewidth=0.8, color="red") + 
  geom_line(data = data.frame(x = c(R0_main_Cpal_peaks_0025, R0_main_Cpal_peaks_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cpal_Tmin_0025, R0_main_Cpal_Tmin_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cpal_Tmax_0025, R0_main_Cpal_Tmax_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_point(data = data.frame(x = c(R0_main_Cpal_Tmin_mean, R0_main_Cpal_peaks_mean, R0_main_Cpal_Tmax_mean), y = 0), aes(x = x, y = y), size=2) + 
  geom_point(data = data.frame(x = c(R0_main_Cpal_peaks_mean), y = 0), aes(x = x, y = y), size=2) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(5,40)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. pipiens pallens")))) +
  theme(axis.text = element_text(size = 14),  
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Cx. restuans

# calculate R0 samples from the trait samples
R0_main_Cres = R0_main_f(biting_fit_pop, infprob_fit_pop, EIP_fit_pop, lf_fit_Cres, omega, sur_fit_Cres, ER, egg_viability_fit_pop, dev_fit_Cres)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Cres <= 0) == ncol(R0_main_Cres))){
  print("lol")
  R0_main_Cres <- R0_main_Cres[-which(rowSums(R0_main_Cres <= 0) == ncol(R0_main_Cres)),]
}

# calculate statistics of relative R0 at each temperature
R0_main_Cres_mean = apply(R0_main_Cres, 2, mean)
R0_main_Cres_0025 = apply(R0_main_Cres, 2, quantile, probs=c(0.025))
R0_main_Cres_0975 = apply(R0_main_Cres, 2, quantile, probs=c(0.975))

# collect information in dataframe and normalize by value at optimal temperature
df_Cres <- data.frame(x = temp, mean = c(R0_main_Cres_mean/max(R0_main_Cres_mean)),
                      lowerCI = c(R0_main_Cres_0025/max(R0_main_Cres_mean)),
                      upperCI = c(R0_main_Cres_0975/max(R0_main_Cres_mean)),
                      Species = c(rep("res.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_main_Cres_peaks = sapply(apply(R0_main_Cres,1,which.max), FUN=function(x)temp[x])
R0_main_Cres_peaks_mean = mean(R0_main_Cres_peaks)
R0_main_Cres_peaks_median =median(R0_main_Cres_peaks)
R0_main_Cres_peaks_0025 = quantile(R0_main_Cres_peaks, probs=c(0.025))
R0_main_Cres_peaks_0975 = quantile(R0_main_Cres_peaks, probs=c(0.975))

# check if any sample has a positive value at lowest temperature which would make 
# Tmin calculations nonsensical
any(R0_main_Cres[,1]>0)

# calculate Tmin for each sample and calculate statistics
R0_main_Cres_Tmin = apply(R0_main_Cres,1,FUN=Tmin)
R0_main_Cres_Tmin_mean = mean(R0_main_Cres_Tmin)
R0_main_Cres_Tmin_median =median(R0_main_Cres_Tmin)
R0_main_Cres_Tmin_0025 = quantile(R0_main_Cres_Tmin, probs=c(0.025))
R0_main_Cres_Tmin_0975 = quantile(R0_main_Cres_Tmin, probs=c(0.975))

# check if any sample has a positive value at highest temperature which would make 
# Tmax calculations nonsensical
any(R0_main_Cres[,451]>0)

# calculate Tmax for each sample and calculate statistics
R0_main_Cres_Tmax = apply(R0_main_Cres,1,FUN=Tmax)
R0_main_Cres_Tmax_mean = mean(R0_main_Cres_Tmax)
R0_main_Cres_Tmax_median =median(R0_main_Cres_Tmax)
R0_main_Cres_Tmax_0025 = quantile(R0_main_Cres_Tmax, probs=c(0.025))
R0_main_Cres_Tmax_0975 = quantile(R0_main_Cres_Tmax, probs=c(0.975))

# plot mean relative R0 temperature response with mean and 95% CIs for 
# Tmin, Topt, Tmax
plot_Cres <- ggplot() +
  geom_line(df_Cres, mapping = aes(x = x, y = mean),linewidth=0.8, color="red") + 
  geom_line(data = data.frame(x = c(R0_main_Cres_peaks_0025, R0_main_Cres_peaks_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cres_Tmin_0025, R0_main_Cres_Tmin_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Cres_Tmax_0025, R0_main_Cres_Tmax_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_point(data = data.frame(x = c(R0_main_Cres_Tmin_mean, R0_main_Cres_peaks_mean, R0_main_Cres_Tmax_mean), y = 0), aes(x = x, y = y), size=2) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(4.9,40)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. restuans")))) +
  theme(axis.text = element_text(size = 14),  
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

# Cx. tarsalis

# calculate R0 samples from the trait samples
R0_main_Ctar = R0_main_f(biting_fit_Ctar, infprob_fit_pop, EIP_fit_Ctar, lf_fit_Ctar, omega, sur_fit_Ctar, ER, egg_viability_fit_pop, dev_fit_Ctar)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Tmin, Tmax, Topt calculations nonsensical
if(any(rowSums(R0_main_Ctar <= 0) == ncol(R0_main_Ctar))){
  print("lol")
  R0_main_Ctar <- R0_main_Ctar[-which(rowSums(R0_main_Ctar <= 0) == ncol(R0_main_Ctar)),]
}

# calculate statistics of relative R0 at each temperature
R0_main_Ctar_mean = apply(R0_main_Ctar, 2, mean)
R0_main_Ctar_0025 = apply(R0_main_Ctar, 2, quantile, probs=c(0.025))
R0_main_Ctar_0975 = apply(R0_main_Ctar, 2, quantile, probs=c(0.975))

# collect information in dataframe and normalize by value at optimal temperature
df_Ctar <- data.frame(x = temp, mean = c(R0_main_Ctar_mean/max(R0_main_Ctar_mean)),
                      lowerCI = c(R0_main_Ctar_0025/max(R0_main_Ctar_mean)),
                      upperCI = c(R0_main_Ctar_0975/max(R0_main_Ctar_mean)),
                      Species = c(rep("tar.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_main_Ctar_peaks = sapply(apply(R0_main_Ctar,1,which.max), FUN=function(x)temp[x])
R0_main_Ctar_peaks_mean = mean(R0_main_Ctar_peaks)
R0_main_Ctar_peaks_median =median(R0_main_Ctar_peaks)
R0_main_Ctar_peaks_0025 = quantile(R0_main_Ctar_peaks, probs=c(0.025))
R0_main_Ctar_peaks_0975 = quantile(R0_main_Ctar_peaks, probs=c(0.975))

# check if any sample has a positive value at lowest temperature which would make 
# Tmin calculations nonsensical
any(R0_main_Ctar[,1]>0)

# calculate Tmin for each sample and calculate statistics
R0_main_Ctar_Tmin = apply(R0_main_Ctar,1,FUN=Tmin)
R0_main_Ctar_Tmin_mean = mean(R0_main_Ctar_Tmin)
R0_main_Ctar_Tmin_median =median(R0_main_Ctar_Tmin)
R0_main_Ctar_Tmin_0025 = quantile(R0_main_Ctar_Tmin, probs=c(0.025))
R0_main_Ctar_Tmin_0975 = quantile(R0_main_Ctar_Tmin, probs=c(0.975))

# check if any sample has a positive value at highest temperature which would make 
# Tmax calculations nonsensical
any(R0_main_Ctar[,451]>0)

# calculate Tmax for each sample and calculate statistics
R0_main_Ctar_Tmax = apply(R0_main_Ctar,1,FUN=Tmax)
R0_main_Ctar_Tmax_mean = mean(R0_main_Ctar_Tmax)
R0_main_Ctar_Tmax_median =median(R0_main_Ctar_Tmax)
R0_main_Ctar_Tmax_0025 = quantile(R0_main_Ctar_Tmax, probs=c(0.025))
R0_main_Ctar_Tmax_0975 = quantile(R0_main_Ctar_Tmax, probs=c(0.975))

# plot mean relative R0 temperature response with mean and 95% CIs for 
# Tmin, Topt, Tmax
plot_Ctar <- ggplot() +
  geom_line(df_Ctar, mapping = aes(x = x, y = mean),linewidth=0.8, color="red") +
  geom_line(data = data.frame(x = c(R0_main_Ctar_peaks_0025, R0_main_Ctar_peaks_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Ctar_Tmin_0025, R0_main_Ctar_Tmin_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_line(data = data.frame(x = c(R0_main_Ctar_Tmax_0025, R0_main_Ctar_Tmax_0975), y = 0), aes(x = x, y = y), linewidth=0.8,alpha=1) + 
  geom_point(data = data.frame(x = c(R0_main_Ctar_Tmin_mean, R0_main_Ctar_peaks_mean, R0_main_Ctar_Tmax_mean), y = 0), aes(x = x, y = y), size=2) + 
  scale_x_continuous(breaks = seq(0,45,5), limits = c(5,40)) +
  theme_bw() +
  ggtitle(expression(paste(italic("Cx. tarsalis")))) +
  theme(axis.text = element_text(size = 14),  
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Combine the relative R0 plots for the six Culex mosquitoes (Figure 3)

plot_list = list(plot_Cpip, plot_Cqui, plot_Cmol, plot_Cpal, plot_Cres, plot_Ctar)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B', 'C', 'D', 'E', 'F'))

y.grob <- textGrob(expression("Transmission suitability " * R[0]^rel), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/R0.pdf", width=8.27, height=9.27)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()

### Comparison to relative R0 models using alternative mosquito abundance approximations

# Alternative model 1 (Eq. S6.2)
M_alt1_f <- function(omega, surJ, ER, a, EV, lf, devJ){
  EFD = ER * a
  surEJ = surJ * EV
  muL = ifelse(surJ>0 & devJ>0, 
               ifelse(surJ<0.99, pmin(devJ * (1-surJ)/(surJ), 1), 
                      pmin(devJ * (1-0.99)/(0.99), 1)), 
               1)
  ifelse(1/(lf * omega * EFD * surEJ + cc) <1,
         omega^2 * EFD * EV * (devJ^2/muL) * lf^2 * (1 - 1/(lf * omega * EFD * surEJ + cc)),
         0)
}

R0_alt1_f = function(a, b, EIP, lf, omega, surJ, ER, EV, devJ){
  muM = 1/(lf+cc)
  (a^2 * b * exp(-EIP*muM) * M_alt1_f(omega, surJ, ER, a, EV, lf, devJ))/(muM)
}

# Alternative model 2 (Eq. S6.4)
M_alt2_f <- function(omega, surJ, ER, a, EV, lf, devJ){
  EFD = ER * a
  surEJ = surJ * EV
  muM = 1/(lf+cc)
  ifelse(1/(lf * omega * EFD * surEJ + cc) <1,
         (omega*devJ/(muM+cc)) * (1 - 1/(lf * omega * EFD * surEJ + cc)),
         0)
}

R0_alt2_f = function(a, b, EIP, lf, omega, surJ, ER, EV, devJ){
  muM = 1/(lf+cc)
  (a^2 * b * exp(-EIP*muM) * M_alt2_f(omega, surJ, ER, a, EV, lf, devJ))/(muM)
}

# Alternative model 3 (Eq. S6.6)
M_alt3_f <- function(omega, surJ, ER, a, EV, lf, devJ){
  EFD = ER * a
  surEJ = surJ * EV
  muM = 1/(lf+cc)
  ifelse(1/(lf * omega * EFD * surEJ + cc) <1,
         (1 - 1/(lf * omega * EFD * surEJ + cc)),
         0)
}

R0_alt3_f = function(a, b, EIP, lf, omega, surJ, ER, EV, devJ){
  muM = 1/(lf+cc)
  (a^2 * b * exp(-EIP*muM) * M_alt3_f(omega, surJ, ER, a, EV, lf, devJ))/(muM)
}

# Alternative model 4 (Eq. S6.7)
M_alt4_f <- function(omega, surJ, ER, a, EV, lf, devJ){
  EFD = ER * a
  EFD * EV * surJ * devJ * lf^2
}

R0_alt4_f = function(a, b, EIP, lf, omega, surJ, ER, EV, devJ){
  muM = 1/(lf+cc)
  (a^2 * b * exp(-EIP*muM) * M_alt4_f(omega, surJ, ER, a, EV, lf, devJ))/(muM)
}

## relative R0 calculations for each model at the example of Culex pipiens
## including calculation of statistics for Topt + plotting

## Alternative model 1

# calculate R0 samples from the trait samples
R0_alt1_Cpip = R0_alt1_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, lf_fit_Cpip, omega, sur_fit_Cpip, ER, egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Topt calculations nonsensical
if(any(rowSums(R0_alt1_Cpip <= 0) == ncol(R0_alt1_Cpip))){
  print("only zeros")
  R0_alt1_Cpip <- R0_alt1_Cpip[-which(rowSums(R0_alt1_Cpip <= 0) == ncol(R0_alt1_Cpip)),]  
}

# calculate statistics of relative R0 at each temperature
R0_alt1_Cpip_mean = apply(R0_alt1_Cpip, 2, mean)
R0_alt1_Cpip_0025 = apply(R0_alt1_Cpip, 2, quantile, probs=c(0.025))
R0_alt1_Cpip_0975 = apply(R0_alt1_Cpip, 2, quantile, probs=c(0.975))

# collect information in dataframe and normalize by value at optimal temperature
df_Cpip_alt1 <- data.frame(x = temp, mean = c(R0_alt1_Cpip_mean/max(R0_alt1_Cpip_mean)),
                      lowerCI = c(R0_alt1_Cpip_0025/max(R0_alt1_Cpip_mean)),
                      upperCI = c(R0_alt1_Cpip_0975/max(R0_alt1_Cpip_mean)),
                      Species = c(rep("pip.",length(temp))))


# calculate Topt for each sample and calculate statistics
R0_alt1_Cpip_peaks = sapply(apply(R0_alt1_Cpip,1,which.max), FUN=function(x)temp[x])
R0_alt1_Cpip_peaks_mean = mean(R0_alt1_Cpip_peaks)
R0_alt1_Cpip_peaks_median =median(R0_alt1_Cpip_peaks)
R0_alt1_Cpip_peaks_0025 = quantile(R0_alt1_Cpip_peaks, probs=c(0.025))
R0_alt1_Cpip_peaks_0975 = quantile(R0_alt1_Cpip_peaks, probs=c(0.975))

## Alternative model 2

# calculate R0 samples from the trait samples
R0_alt2_Cpip = R0_alt2_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, lf_fit_Cpip, omega, sur_fit_Cpip, ER, egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Topt calculations nonsensical
if(any(rowSums(R0_alt2_Cpip <= 0) == ncol(R0_alt2_Cpip))){
  print("only zeros")
  R0_alt2_Cpip <- R0_alt2_Cpip[-which(rowSums(R0_alt2_Cpip <= 0) == ncol(R0_alt2_Cpip)),]  
}

# calculate statistics of relative R0 at each temperature
R0_alt2_Cpip_mean = apply(R0_alt2_Cpip, 2, mean)
R0_alt2_Cpip_0025 = apply(R0_alt2_Cpip, 2, quantile, probs=c(0.025))
R0_alt2_Cpip_0975 = apply(R0_alt2_Cpip, 2, quantile, probs=c(0.975))

# collect information in dataframe and normalize by value at optimal temperature
df_Cpip_alt2 <- data.frame(x = temp, mean = c(R0_alt2_Cpip_mean/max(R0_alt2_Cpip_mean)),
                      lowerCI = c(R0_alt2_Cpip_0025/max(R0_alt2_Cpip_mean)),
                      upperCI = c(R0_alt2_Cpip_0975/max(R0_alt2_Cpip_mean)),
                      Species = c(rep("pip.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_alt2_Cpip_peaks = sapply(apply(R0_alt2_Cpip,1,which.max), FUN=function(x)temp[x])
R0_alt2_Cpip_peaks_mean = mean(R0_alt2_Cpip_peaks)
R0_alt2_Cpip_peaks_median =median(R0_alt2_Cpip_peaks)
R0_alt2_Cpip_peaks_0025 = quantile(R0_alt2_Cpip_peaks, probs=c(0.025))
R0_alt2_Cpip_peaks_0975 = quantile(R0_alt2_Cpip_peaks, probs=c(0.975))

## Alternative model 3

# calculate R0 samples from the trait samples
R0_alt3_Cpip = R0_alt3_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, lf_fit_Cpip, omega, sur_fit_Cpip, ER, egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Topt calculations nonsensical
if(any(rowSums(R0_alt3_Cpip <= 0) == ncol(R0_alt3_Cpip))){
  print("only zeros")
  R0_alt3_Cpip <- R0_alt3_Cpip[-which(rowSums(R0_alt3_Cpip <= 0) == ncol(R0_alt3_Cpip)),]  
}

# calculate statistics of relative R0 at each temperature
R0_alt3_Cpip_mean = apply(R0_alt3_Cpip, 2, mean)
R0_alt3_Cpip_0025 = apply(R0_alt3_Cpip, 2, quantile, probs=c(0.025))
R0_alt3_Cpip_0975 = apply(R0_alt3_Cpip, 2, quantile, probs=c(0.975))

# collect information in dataframe and normalize by value at optimal temperature
df_Cpip_alt3 <- data.frame(x = temp, mean = c(R0_alt3_Cpip_mean/max(R0_alt3_Cpip_mean)),
                      lowerCI = c(R0_alt3_Cpip_0025/max(R0_alt3_Cpip_mean)),
                      upperCI = c(R0_alt3_Cpip_0975/max(R0_alt3_Cpip_mean)),
                      Species = c(rep("pip.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_alt3_Cpip_peaks = sapply(apply(R0_alt3_Cpip,1,which.max), FUN=function(x)temp[x])
R0_alt3_Cpip_peaks_mean = mean(R0_alt3_Cpip_peaks)
R0_alt3_Cpip_peaks_median =median(R0_alt3_Cpip_peaks)
R0_alt3_Cpip_peaks_0025 = quantile(R0_alt3_Cpip_peaks, probs=c(0.025))
R0_alt3_Cpip_peaks_0975 = quantile(R0_alt3_Cpip_peaks, probs=c(0.975))

## Alternative model 4

# calculate R0 samples from the trait samples
R0_alt4_Cpip = R0_alt4_f(biting_fit_Cpip, infprob_fit_pop, EIP_fit_Cpip, lf_fit_Cpip, omega, sur_fit_Cpip, ER, egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Topt calculations nonsensical
if(any(rowSums(R0_alt4_Cpip <= 0) == ncol(R0_alt4_Cpip))){
  print("only zeros")
  R0_alt4_Cpip <- R0_alt4_Cpip[-which(rowSums(R0_alt4_Cpip <= 0) == ncol(R0_alt4_Cpip)),]  
}

# calculate statistics of relative R0 at each temperature
R0_alt4_Cpip_mean = apply(R0_alt4_Cpip, 2, mean)
R0_alt4_Cpip_0025 = apply(R0_alt4_Cpip, 2, quantile, probs=c(0.025))
R0_alt4_Cpip_0975 = apply(R0_alt4_Cpip, 2, quantile, probs=c(0.975))

# collect information in dataframe and normalize by value at optimal temperature
df_Cpip_alt4 <- data.frame(x = temp, mean = c(R0_alt4_Cpip_mean/max(R0_alt4_Cpip_mean)),
                      lowerCI = c(R0_alt4_Cpip_0025/max(R0_alt4_Cpip_mean)),
                      upperCI = c(R0_alt4_Cpip_0975/max(R0_alt4_Cpip_mean)),
                      Species = c(rep("pip.",length(temp))))

# calculate Topt for each sample and calculate statistics
R0_alt4_Cpip_peaks = sapply(apply(R0_alt4_Cpip,1,which.max), FUN=function(x)temp[x])
R0_alt4_Cpip_peaks_mean = mean(R0_alt4_Cpip_peaks)
R0_alt4_Cpip_peaks_median =median(R0_alt4_Cpip_peaks)
R0_alt4_Cpip_peaks_0025 = quantile(R0_alt4_Cpip_peaks, probs=c(0.025))
R0_alt4_Cpip_peaks_0975 = quantile(R0_alt4_Cpip_peaks, probs=c(0.975))

## Collect mean temperature response of the different relative R0 models in one
## dataframe
df_Cpip_alt = data.frame(x = rep(temp, 5), 
                         model = c(rep("main",length(temp)),
                                   rep("model 1",length(temp)),
                                   rep("model 2",length(temp)),
                                   rep("model 3",length(temp)),
                                   rep("model 4",length(temp))),
                         mean = c(R0_main_Cpip_mean/max(R0_main_Cpip_mean),
                                  R0_alt1_Cpip_mean/max(R0_alt1_Cpip_mean),
                                  R0_alt2_Cpip_mean/max(R0_alt2_Cpip_mean),
                                  R0_alt3_Cpip_mean/max(R0_alt3_Cpip_mean),
                                  R0_alt4_Cpip_mean/max(R0_alt4_Cpip_mean)))

## Plot mean temperature response of the different relative R0 models
plot_Cpip_alt <- ggplot(df_Cpip_alt, aes(x = x, y = mean, color=model)) +
  geom_line(linewidth = 0.8) + 

  scale_x_continuous(breaks = seq(0,45,5), limits = c(5,40)) +
  theme_bw() +
  labs(title = expression(paste(italic("Cx. pipiens"))),
       x = "Temperature (°C)",
       y = expression("Transmission suitability " * R[0]^rel), 
       color = "Mosquito model") + 
  scale_color_manual(values = c("red", "orange", "black", "blue", "#785EF0")) +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 10),
        legend.position = c(.2,.7),
        legend.text = element_text(size = 10),
        plot.title = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Collect statistics about Topt for the different model in one dataframe
df_alt_peaks = data.frame(label = c("5main","4model1","3model2","2model3","1model4"), 
                          mean = c(R0_main_Cpip_peaks_mean, R0_alt1_Cpip_peaks_mean, 
                                   R0_alt2_Cpip_peaks_mean, R0_alt3_Cpip_peaks_mean, 
                                   R0_alt4_Cpip_peaks_mean), 
                          lower = c(R0_main_Cpip_peaks_0025, R0_alt1_Cpip_peaks_0025,
                                    R0_alt2_Cpip_peaks_0025, R0_alt3_Cpip_peaks_0025,
                                    R0_alt4_Cpip_peaks_0025), 
                          upper = c(R0_main_Cpip_peaks_0975, R0_alt1_Cpip_peaks_0975,
                                    R0_alt2_Cpip_peaks_0975, R0_alt3_Cpip_peaks_0975,
                                    R0_alt4_Cpip_peaks_0975))

## Plot statistics about Topt for the different model
plot_Cpip_alt_peaks <- ggplot(df_alt_peaks, aes(x=label, y=mean, ymin=lower, ymax=upper, color=label)) +
  geom_pointrange(linewidth=0.8) +
  coord_flip() +
  guides(color = "none") +
  theme_bw() +
  scale_color_manual(values = c("#785EF0", "blue", "black", "orange", "red")) +
  scale_x_discrete(labels = c("model4","model3","model2","model1","main")) +
  labs(title = expression(paste("optimal Temperature for ", italic("Cx. pip. ") * R[0]^rel)),
       y = "Temperature (°C)") + 
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), # Adjust axis title size
        axis.title.y = element_blank(),
        #legend.key.size = unit(1.5, "lines"),   # Larger legend key size
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Combine relative R0 and Topt plots for different models (Figure 4)

plot_list = list(plot_Cpip_alt, plot_Cpip_alt_peaks)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B'))

#pdf("Figures/R0_compare_models_Cpip.pdf", width=8.3, height=4)
plot_grid
#dev.off()

### relative R0 calculations including Topt with different levels of 
### between-experiment variability for mosquito infection probability and EIP 
### at the example of Culex pipiens (SI3)

## increased between-experiment variability 

# calculate R0 samples from the trait samples
R0_main_Cpip_inc = R0_main_f(biting_fit_Cpip, infprob_fit_pop_inc, EIP_fit_Cpip_inc, lf_fit_Cpip, omega, sur_fit_Cpip, ER, egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_inc <= 0) == ncol(R0_main_Cpip_inc))){
  print("only zeros")
  R0_main_Cpip_inc <- R0_main_Cpip_inc[-which(rowSums(R0_main_Cpip_inc <= 0) == ncol(R0_main_Cpip_inc)),]  
}

# calculate statistics of relative R0 at each temperature
R0_main_Cpip_inc_mean = apply(R0_main_Cpip_inc, 2, mean)
R0_main_Cpip_inc_0025 = apply(R0_main_Cpip_inc, 2, quantile, probs=c(0.025))
R0_main_Cpip_inc_0975 = apply(R0_main_Cpip_inc, 2, quantile, probs=c(0.975))

# calculate Topt for each sample and calculate statistics
R0_main_Cpip_inc_peaks = sapply(apply(R0_main_Cpip_inc,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_inc_peaks_mean = mean(R0_main_Cpip_inc_peaks)
R0_main_Cpip_inc_peaks_median = median(R0_main_Cpip_inc_peaks)
R0_main_Cpip_inc_peaks_0025 = quantile(R0_main_Cpip_inc_peaks, probs=c(0.025))
R0_main_Cpip_inc_peaks_0975 = quantile(R0_main_Cpip_inc_peaks, probs=c(0.975))

## reduced between-experiment variability 

# calculate R0 samples from the trait samples
R0_main_Cpip_red = R0_main_f(biting_fit_Cpip, infprob_fit_pop_red, EIP_fit_Cpip_red, lf_fit_Cpip, omega, sur_fit_Cpip, ER, egg_viability_fit_pop, dev_fit_Cpip)

# check for and remove any sample relative R0 samples has only zeros across whole 
# temperature range which would make Topt calculations nonsensical
if(any(rowSums(R0_main_Cpip_red <= 0) == ncol(R0_main_Cpip_red))){
  print("only zeros")
  R0_main_Cpip_red <- R0_main_Cpip_red[-which(rowSums(R0_main_Cpip_red <= 0) == ncol(R0_main_Cpip_red)),]  
}

# calculate statistics of relative R0 at each temperature
R0_main_Cpip_red_mean = apply(R0_main_Cpip_red, 2, mean)
R0_main_Cpip_red_0025 = apply(R0_main_Cpip_red, 2, quantile, probs=c(0.025))
R0_main_Cpip_red_0975 = apply(R0_main_Cpip_red, 2, quantile, probs=c(0.975))

# calculate Topt for each sample and calculate statistics
R0_main_Cpip_red_peaks = sapply(apply(R0_main_Cpip_red,1,which.max), FUN=function(x)temp[x])
R0_main_Cpip_red_peaks_mean = mean(R0_main_Cpip_red_peaks)
R0_main_Cpip_red_peaks_median = median(R0_main_Cpip_red_peaks)
R0_main_Cpip_red_peaks_0025 = quantile(R0_main_Cpip_red_peaks, probs=c(0.025))
R0_main_Cpip_red_peaks_0975 = quantile(R0_main_Cpip_red_peaks, probs=c(0.975))

# collect Topt statistics for the different between-experiment variability 
# scenarios in one dataframe
df_Cpip_variances_peaks = data.frame(label = c("3reduced","2main","1increased"), 
                          mean = c(R0_main_Cpip_red_peaks_mean, R0_main_Cpip_peaks_mean, 
                                   R0_main_Cpip_inc_peaks_mean), 
                          lower = c(R0_main_Cpip_red_peaks_0025, R0_main_Cpip_peaks_0025,
                                    R0_main_Cpip_inc_peaks_0025), 
                          upper = c(R0_main_Cpip_red_peaks_0975, R0_main_Cpip_peaks_0975,
                                    R0_main_Cpip_inc_peaks_0975))

# Plot Topt statistics for the different between-experiment variability (Figure SI4.5)
plot_Cpip_variances_peaks <- ggplot(df_Cpip_variances_peaks, aes(x=label, y=mean, ymin=lower, ymax=upper, color=label)) +
  geom_pointrange(linewidth=0.8) +
  coord_flip() +
  guides(color = "none") +
  theme_bw() +
  scale_color_manual(values = c("red","orange", "black")) +
  scale_x_discrete(labels = c("small var.","moderate (main) var.","large var.")) +
  labs(title = expression(paste("optimal Temperature for ", italic("Cx. pipiens ") * R[0]^rel)),
       y = "Temperature (°C)") + 
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 11),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        plot.title = element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#pdf("Figures/R0_compare_variabilities_Cpip.pdf", width=4.5, height=3)
plot_Cpip_variances_peaks 
#dev.off()

