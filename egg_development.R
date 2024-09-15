### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to fit temperature-dependent functions for mosquito egg
### development rates (from egg laying to egg hatch) using a multi-species 
### Bayesian hierarchical model implemented in STAN 
###
### The script also includes the necessary code to reproduce Figure 6
###

## load libraries
library(rstan)
library(tidyverse)
library(mcmcplots)
library(bayesplot)
library(posterior)
library(cowplot)
library(grid)
library(gridExtra)

# enable parallel computation and allow STAN to automatically overwrite models
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE) 

# set working directory
## setwd("...")

## load data
data_egg_dev <- read.csv("data/egg_development_data.csv",sep=",")

# the original data is given as development times --> transform to rate=1/time
data_egg_dev$trait <- 1/data_egg_dev$trait

# create a separate dataframe for each species
data_egg_dev_Cpip <- data_egg_dev %>% filter(species_id=="Cpip")
data_egg_dev_Cqui <- data_egg_dev %>% filter(species_id=="Cqui")
data_egg_dev_Cres <- data_egg_dev %>% filter(species_id=="Cres") 
data_egg_dev_Cmel <- data_egg_dev %>% filter(species_id=="Cmel")
data_egg_dev_Cmol <- data_egg_dev %>% filter(species_id=="Cmol")
data_egg_dev_Cpal <- data_egg_dev %>% filter(species_id=="Cpal")
data_egg_dev_Cthei <- data_egg_dev %>% filter(species_id=="Cthei")

## prepare input data for STAN model

# temp. points on which to calculate generated quantities
steps = 0.1
temp <- seq(0,45,steps) 
# a vector storing integer ids for species each observation belongs to
species_id = as.integer(as.factor(data_egg_dev$species_id)) 
# a vector storing integer ids for experiment each observation belongs to
exp_id = as.integer(as.factor(data_egg_dev$experiment_id)) 

# data list collecting all information needed to initialize the STAN model
data <- list(Nspecies = length(unique(data_egg_dev$species_id)), # number of species
             Nexp = length(unique(data_egg_dev$experiment_id)), # number of experiments
             Nobs = nrow(data_egg_dev), # total number of observations
             x = data_egg_dev$temperature, # predictor (temperature)
             y = data_egg_dev$trait, # outcome (development rate)
             species_id = species_id, # species id of each observation
             exp_id = exp_id, # experiment id of each observation
             N_new = length(temp), # number of temp. points for generated quantities
             x_new = temp) # temp. points for generated quantities

## Fitting the STAN model or reading a fitted model

# fitting
#fit <- stan(file = 'stan_models/egg_development_model.stan',
#            data = data, iter=4000, chains=4,
#            control = list(adapt_delta=0.99, max_treedepth=12))

# save fitted model 
#saveRDS(fit, "model_fits/egg_development_fit.rds")

# read fitted model
fit <- readRDS("model_fits/egg_development_fit.rds")

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

### Visualizations

## Basic plots of mean model fit (+ 95% credible interval) for each species

# Extract model fit summary (generated quantity)
f_new <- summary(fit, pars="f_new")$summary # model fit summary for each temperature (n=451) and each species (n=7) --> 3157 rows
f_new_mean <- matrix(f_new[,1], nrow=data$N_new) # mean of model fits
f_new_0025 <- matrix(f_new[,4], nrow=data$N_new) # 2.5% quantile of model fits
f_new_0975 <- matrix(f_new[,8], nrow=data$N_new) # 97.5% quantile of model fits

# Function to make basic plots of mean model fit (+ 95% credible interval) for each species
plot.species <- function(i, x, y, title){
  plot(temp,f_new_mean[,i],type="l",ylim=c(0,1.5),
       xlab="Temperature (°C)",
       ylab="Egg development rate",
       lwd=2,cex.axis=1,cex.lab=1, main=title)
  lines(temp,f_new_0025[,i],lty=2,lwd=2)
  lines(temp,f_new_0975[,i],lty=2,lwd=2)
  points(x,y,cex=1,pch=19)
} 

# plot estimates of each species' expected response
#Cmel plot
plot.species(1, data_egg_dev_Cmel$temperature,data_egg_dev_Cmel$trait, title="Culiseta melanura")

#Cmol plot
plot.species(2, data_egg_dev_Cmol$temperature,data_egg_dev_Cmol$trait, title="Culex pip. molestus")

#Cpal plot
plot.species(3, data_egg_dev_Cpal$temperature,data_egg_dev_Cpal$trait, title="Culex pip. pallens")

#Cpip plot
plot.species(4, data_egg_dev_Cpip$temperature,data_egg_dev_Cpip$trait, title="Culex pipiens")

#Cqui plot
plot.species(5, data_egg_dev_Cqui$temperature,data_egg_dev_Cqui$trait, title="Culex quinque")

#Cres plot
plot.species(6, data_egg_dev_Cres$temperature,data_egg_dev_Cres$trait, title="Culex restuans")

#Cthei plot
plot.species(7, data_egg_dev_Cthei$temperature,data_egg_dev_Cthei$trait, title="Culex theileri")

## Visualization of expected temperature response of a new species generated by sampling from hierachical prior
f_new_spec <- summary(fit, pars="f_new_spec")$summary # extract summary for expected response of a new species
plot(temp,f_new_spec[,1],type="l",ylim=c(0,1.5),ylab="Egg development rate",
     lwd=2,cex.axis=1,cex.lab=1, main="New species")
lines(temp,f_new_spec[,4],lty=2,lwd=2)
lines(temp,f_new_spec[,8],lty=2,lwd=2)

### more sophisticated plots for article

## Plots of mean model fit (+ 95% credible interval) for each species

#Cmel plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,1], 
                     lowerCI = f_new_0025[,1], upperCI = f_new_0975[,1])
plot1 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_egg_dev_Cmel, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culiseta melanura"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5))

#Cmol plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,2], 
                     lowerCI = f_new_0025[,2], upperCI = f_new_0975[,2])
plot2 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_egg_dev_Cmol, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens molestus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5))

#Cpal plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,3], 
                     lowerCI = f_new_0025[,3], upperCI = f_new_0975[,3])
plot3 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_egg_dev_Cpal, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens pallens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5))

#Cpip plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,4], 
                     lowerCI = f_new_0025[,4], upperCI = f_new_0975[,4])
plot4 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_egg_dev_Cpip, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5))

#Cqui plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,5], 
                     lowerCI = f_new_0025[,5], upperCI = f_new_0975[,5])
plot5 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_egg_dev_Cqui, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex quinquefasciatus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5))

#Cres plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,6], 
                     lowerCI = f_new_0025[,6], upperCI = f_new_0975[,6])
plot6 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_egg_dev_Cres, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex restuans"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5))

#Cthei plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,7], 
                     lowerCI = f_new_0025[,7], upperCI = f_new_0975[,7])
plot7 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_egg_dev_Cthei, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex theileri"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1.5)) 

# Plot six Culex species together (Figure 6)
plot_list = list(plot4, plot5, plot2, plot3, plot6)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B', 'C', 'D', 'E', 'F'))

y.grob <- textGrob(expression(paste("Egg to larva development rate (days"^{-1},")")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/egg_development.pdf", width=8.27, height=9.27)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()

