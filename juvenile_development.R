### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to fit temperature-dependent functions for mosquito juvenile
### development rates (from egg hatch to adult emergence) using a multi-species 
### Bayesian hierarchical model implemented in STAN 
###
### The script also includes the necessary code to reproduce Figures 5, SI5.1,
### SI5.2
###

## load libraries
library(rstan)
library(tidyverse)
library(mcmcplots)
library(bayesplot)
library(posterior)
library(grid)
library(gridExtra)
library(MASS)

# enable parallel computation and allow STAN to automatically overwrite models
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE) 

# set working directory
## setwd("...")

## load data
data_dev <- read.csv("data/juvenile_development_data.csv",sep=",")

# the original data is given as development times --> transform to rate=1/time
data_dev$trait <- 1/data_dev$trait

# create a separate dataframe for each species
data_dev_Asol <- data_dev %>% filter(species_id=="Asol")
data_dev_Atri <- data_dev %>% filter(species_id=="Atri")
data_dev_Cino <- data_dev %>% filter(species_id=="Cino")
data_dev_Cmel <- data_dev %>% filter(species_id=="Cmel")
data_dev_Cmol <- data_dev %>% filter(species_id=="Cmol")
data_dev_Cpal <- data_dev %>% filter(species_id=="Cpal")
data_dev_Cpip <- data_dev %>% filter(species_id=="Cpip")
data_dev_Cqui <- data_dev %>% filter(species_id=="Cqui")
data_dev_Cres <- data_dev %>% filter(species_id=="Cres") 
data_dev_Csal <- data_dev %>% filter(species_id=="Csal")
data_dev_Ctar <- data_dev %>% filter(species_id=="Ctar")

## prepare input data for STAN model

# temp. points on which to calculate generated quantities
steps = 0.1
temp <- seq(0,45,steps) 
# a vector storing integer ids for species each observation belongs to
species_id = as.integer(as.factor(data_dev$species_id)) 
# a vector storing integer ids for experiment each observation belongs to
exp_id = as.integer(as.factor(data_dev$experiment_id)) 

# data list collecting all information needed to initialize the STAN model
data <- list(Nspecies = length(unique(data_dev$species_id)), # number of species
             Nexp = length(unique(data_dev$experiment_id)), # number of experiments
             Nobs = nrow(data_dev), # total number of observations
             x = data_dev$temperature, # predictor (temperature)
             y = data_dev$trait, # outcome (development rate)
             species_id = species_id, # species id of each observation
             exp_id = exp_id, # experiment id of each observation
             N_new = length(temp), # number of temp. points for generated quantities
             x_new = temp) # temp. points for generated quantities

## Fitting the STAN model or reading a fitted model

# fitting
#fit <- stan(file = 'stan_models/juvenile_development_model.stan',
#           data = data, iter=4000, warmup=2000, chains=4, 
#           control = list(adapt_delta = 0.99, max_treedepth = 12))

# save fitted model 
#saveRDS(fit, "model_fits/juvenile_development_fit.rds")

# read fitted model
fit <- readRDS("model_fits/juvenile_development_fit.rds")

## inspect the STAN model fit 

# print model summary
print(fit)

# print between species and experiment variability estimates 
print(fit,pars=c("sigma_exp_q","sigma_exp_Tmin","sigma_exp_Tmax"),probs = c(0.025,0.975))
print(fit,pars=c("sigma_q","sigma_Tmin","sigma_Tmax"),probs = c(0.025,0.975))

# print species level parameter estimates 
print(fit,pars=c("q"),probs = c(0.025,0.975))
print(fit,pars=c("Tmax"),probs = c(0.025,0.975))
print(fit,pars=c("Tmin"),probs = c(0.025,0.975))
print(fit,pars=c("Topt"),probs = c(0.025,0.975))

# look at pairs plots to detect potential sampling problems
pairs(fit, pars=c("mu_q","sigma_q","sigma_exp_q",
                  "mu_Tmin","sigma_Tmin","sigma_exp_Tmin",
                  "mu_Tmax","sigma_Tmax","sigma_exp_Tmax")) # population-level mean and between species/experiment varibility
mcmc_pairs(fit, pars=c("sigma_exp_q","sigma_exp_Tmin","sigma_exp_Tmax"), 
           diag_fun="dens", off_diag_fun = "hex") # different visualization using the bayesplot package

# look at density- and traceplots to detect potential sampling problems
s = as.array(fit)
mcmc <- do.call(mcmc.list, plyr:::alply(s[, , -(length(s[1, 1, ]))], 2, as.mcmc))
denplot(mcmc, parms = c("mu_q","mu_Tmin","mu_Tmax","sigma_q","sigma_Tmin","sigma_Tmax")) 
denplot(mcmc, parms = c("sigma_exp_q","sigma_exp_Tmin","sigma_exp_Tmax")) 
traplot(mcmc, parms = c("mu_q","sigma_q","mu_Tmin","sigma_Tmin","mu_Tmax","sigma_Tmax")) 
traplot(mcmc, parms = c("sigma_exp_q","sigma_exp_Tmin","sigma_exp_Tmax")) 

### Fit gamma distributions to the posterior distribution samples for between-experiment
### and between-species standard deviations. These are used as informative priors for 
### the same parameters for other traits (see manuscript)

# extract posterior samples for each parameter and combine samples in dataframe
sigma_Tmax_samples <- rstan::extract(fit, permuted=T)$sigma_Tmax
sigma_exp_Tmax_samples <- rstan::extract(fit, permuted=T)$sigma_exp_Tmax
sigma_Tmin_samples <- rstan::extract(fit, permuted=T)$sigma_Tmin
sigma_exp_Tmin_samples <- rstan::extract(fit, permuted=T)$sigma_exp_Tmin
sigma_q_samples <- rstan::extract(fit, permuted=T)$sigma_q
sigma_exp_q_samples <- rstan::extract(fit, permuted=T)$sigma_exp_q

sigma_params_samples <- data.frame(Tmin = sigma_Tmin_samples,
                                   Tmin_exp = sigma_exp_Tmin_samples,
                                   Tmax = sigma_Tmax_samples,
                                   Tmax_exp = sigma_exp_Tmax_samples,
                                   q = sigma_q_samples,
                                   q_exp = sigma_exp_q_samples)

# fit gamma distributions for each parameter and print output
gamma_fit <- apply(sigma_params_samples, 2, function(df) fitdistr(df, "gamma", lower=c(0,0))$estimate)
gamma_fit

### Visualizations

## Basic plots of mean model fit (+ 95% credible interval) for each species

# Extract model fit summary (generated quantity)
f_new <- summary(fit, pars="f_new")$summary # model fit summary for each temperature (n=451) and each species (n=11) --> 4961 rows
f_new_mean <- matrix(f_new[,1], nrow=data$N_new) # mean of model fits
f_new_0025 <- matrix(f_new[,4], nrow=data$N_new) # 2.5% quantile of model fits
f_new_0975 <- matrix(f_new[,8], nrow=data$N_new) # 97.5% quantile of model fits

# Function to make basic plots of mean model fit (+ 95% credible interval) for each species
plot.species <- function(i, x, y, title){
  plot(temp,f_new_mean[,i],type="l",ylim=c(0,0.25),
       xlab="Temperature (°C)",
       ylab="Larva to Adult development rate",
       lwd=2,cex.axis=1,cex.lab=1, main=title)
  lines(temp,f_new_0025[,i],lty=2,lwd=2)
  lines(temp,f_new_0975[,i],lty=2,lwd=2)
  points(x,y,cex=1,pch=19)
} 

# plot estimates of each species' expected response
#Asol plot
plot.species(1, data_dev_Asol$temperature,data_dev_Asol$trait, title="Aedes sollicitans")

#Atri plot
plot.species(2, data_dev_Atri$temperature,data_dev_Atri$trait, title="Aedes triseriatus")

#Cino plot
plot.species(3, data_dev_Cino$temperature,data_dev_Cino$trait, title="Culiseta inornata")

#Cmel plot
plot.species(4, data_dev_Cmel$temperature,data_dev_Cmel$trait, title="Culiseta melanura")

#Cmol plot
plot.species(5, data_dev_Cmol$temperature,data_dev_Cmol$trait, title="Culex pip. molestus")

#Cpal plot
plot.species(6, data_dev_Cpal$temperature,data_dev_Cpal$trait, title="Culex pip. pallens")

#Cpip plot
plot.species(7, data_dev_Cpip$temperature,data_dev_Cpip$trait, title="Culex pipiens")

#Cqui plot
plot.species(8, data_dev_Cqui$temperature,data_dev_Cqui$trait, title="Culex quinque")

#Cres plot
plot.species(9, data_dev_Cres$temperature,data_dev_Cres$trait, title="Culex restuans")

#Csal plot
plot.species(10, data_dev_Csal$temperature,data_dev_Csal$trait, title="Culex salinarius")

#Ctar plot
plot.species(11, data_dev_Ctar$temperature,data_dev_Ctar$trait, title="Culex tarsalis")

## Visualization of expected temperature response of a new species generated by sampling from hierarchical prior
f_new_spec <- summary(fit, pars="f_new_species")$summary # extract summary for expected response of a new species
plot(temp,f_new_spec[,1],type="l",ylim=c(0,0.25),ylab="Larva to Adult development rate",
          lwd=2,cex.axis=1,cex.lab=1, main="New species")
lines(temp,f_new_spec[,4],lty=2,lwd=2)
lines(temp,f_new_spec[,8],lty=2,lwd=2)

### more sophisticated plots for article

## Plots of mean model fit (+ 95% credible interval) for each species

# Asol plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,1], 
                     lowerCI = f_new_0025[,1], upperCI = f_new_0975[,1])
plot1 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Asol, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Aedes sollicitans"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Atri plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,2], 
                     lowerCI = f_new_0025[,2], upperCI = f_new_0975[,2])
plot2 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Atri, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Aedes triseriatus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Cino plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,3], 
                     lowerCI = f_new_0025[,3], upperCI = f_new_0975[,3])
plot3 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cino, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culiseta inornata"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Cmel plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,4], 
                     lowerCI = f_new_0025[,4], upperCI = f_new_0975[,4])
plot4 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cmel, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culiseta melanura"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Cmol plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,5], 
                     lowerCI = f_new_0025[,5], upperCI = f_new_0975[,5])
plot5 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cmol, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens molestus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Cpal plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,6], 
                     lowerCI = f_new_0025[,6], upperCI = f_new_0975[,6])
plot6 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cpal, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens pallens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Cpip plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,7], 
                     lowerCI = f_new_0025[,7], upperCI = f_new_0975[,7])
plot7 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cpip, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Cqui plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,8], 
                     lowerCI = f_new_0025[,8], upperCI = f_new_0975[,8])
plot8 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cqui, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex quinquefasciatus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Cres plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,9], 
                     lowerCI = f_new_0025[,9], upperCI = f_new_0975[,9])
plot9 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cres, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex restuans"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Csal plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,10], 
                     lowerCI = f_new_0025[,10], upperCI = f_new_0975[,10])
plot10 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Csal, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex salinarius"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18))

# Ctar plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,11], 
                     lowerCI = f_new_0025[,11], upperCI = f_new_0975[,11])
plot11 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Ctar, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex tarsalis"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.18)) 

# Plot six Culex species together (Figure 5)
plot_list = list(plot7, plot8, plot5, plot6, plot9, plot11)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                              align = "h", axis = "b", labels = c('A', 'B', 'C', 'D', 'E', 'F'))

y.grob <- textGrob(expression(paste("Larva to adult development rate (days"^{-1},")")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/juvenile_development.pdf", width=8.27, height=9.27)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()

## Visualization of the mean model fit (+ 95% credible interval) for Cx. pipiens along with experiment-specific mean fits (SI4)
f_new_exp <- summary(fit, pars="f_new_exp")$summary # extract summary for experiment-specific fits
f_new_mean_exp <- matrix(f_new_exp[,1],nrow=data$N_new) # mean of experiment-specific fits
# these quantities are only sensible for experiments (exp_ids) conducted on Cx. pipiens

# prepare dataframe for ggplot that contains the mean fits for each experiment conducted on Cx. pipiens
helper = cbind(data.frame(x = temp),f_new_mean_exp) 
helper = helper %>% dplyr::select(1,unique(data_dev_Cpip$experiment_id)+1) %>%
  pivot_longer(cols=2:14,
               names_to = "exp_id",
               values_to = "trait")

# Figure SI5.1
plot_Cpip_studies <- ggplot() + 
  geom_line(data = data.frame(x = temp, y = f_new_mean[,7]), aes(x = x, y = y), color = "black", linewidth = 0.7) +
  geom_ribbon(data = data.frame(x = temp, ymin = f_new_0025[,7], ymax = f_new_0975[,7]), 
              aes(x = x, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.07) +
  geom_line(data = helper, aes(x = x, y = trait, group = exp_id, color = exp_id), linetype = "dashed", linewidth = 0.7) +
  geom_point(data = data_dev_Cpip, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 1.5, shape = 19) +
  labs(x = "Temperature (°C)", y = expression(paste("Larva to adult development rate (days"^{-1},")")), 
       title= expression(paste(italic("Culex pipiens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  coord_cartesian(ylim = c(0, 0.18))

#pdf("Figures/juvenile_development_Cpip_experiments.pdf", width=4.1, height=3.1)
plot_Cpip_studies
#dev.off()

## Comparison of model fits to model fit ignoring experiment identities/variability (SI4)

# Fit or load alternative STAN model that ignores experiment identity effects
#fit2 <- stan(file = 'stan_models/juvenile_development_model_no_experiment_variability.stan',
#            data = data, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99, max_treedepth = 12))
#saveRDS(fit2, "model_fits/juvenile_development_fit_no_experiment_variability.rds.rds")
fit2 <- readRDS("model_fits/juvenile_development_fit_no_experiment_variability.rds")

# Extract model fit2 summary (generated quantity)
f_new2 <- summary(fit2, pars="f_new")$summary # extract model fit2 summary 
f_new_mean2 <- matrix(f_new2[,1],nrow=data$N_new) # mean of model fit
f_new_0025_2 <- matrix(f_new2[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_0975_2 <- matrix(f_new2[,8],nrow=data$N_new) # 97.5% quantile of model fit

## Plots of mean model fit (+ 95% credible interval) for Asol and Cqui for model with and without experiment variability

# Plot for Asol with experiment variability (main model)
temp_df = data.frame(temp = temp, mean = f_new_mean[,1], 
                     lowerCI = f_new_0025[,1], upperCI = f_new_0975[,1])
plot_Asol = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Asol, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Ae. sol. with between-experiment var."))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.25))

# Plot for Asol  without experiment variability 

temp_df = data.frame(temp = temp, mean = f_new_mean2[,1], 
                     lowerCI = f_new_0025_2[,1], upperCI = f_new_0975_2[,1])
plot_Asol_alt = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Asol, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Ae. sol. w/o between-experiment var."))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.25))

# Plot for Cqui with experiment variability (main model)

temp_df = data.frame(temp = temp, mean = f_new_mean[,8], 
                     lowerCI = f_new_0025[,8], upperCI = f_new_0975[,8])
plot_Cqui = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cqui, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Cx. quin. with between-experiment var."))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  coord_cartesian(ylim = c(0, 0.18))

# Plot for Cqui without experiment variability 

temp_df = data.frame(temp = temp, mean = f_new_mean2[,8], 
                     lowerCI = f_new_0025_2[,8], upperCI = f_new_0975_2[,8])
plot_Cqui_alt = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_dev_Cqui, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Cx. quin. w/o between-experiment var."))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  coord_cartesian(ylim = c(0, 0.18))

# Combine model comparison plots together (Figure SI5.2)
plot_list = list(plot_Cqui, plot_Cqui_alt, plot_Asol, plot_Asol_alt)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B', 'C', 'D', 'E', 'F'))

y.grob <- textGrob(expression(paste("Larva to adult development rate (days"^{-1},")")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/juvenile_development_with_without_experiment_variability.pdf", width=8.27, height=6.3)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()

