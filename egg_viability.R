### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to fit temperature-dependent functions for mosquito egg
### viability using a multi-species Bayesian hierarchical model implemented in STAN 
###
### The script also includes the necessary code to reproduce Figure 9
###

## load libraries
library(rstan)
library(tidyverse)
library(mcmcplots)
library(bayesplot)
library(posterior)
library(grid)
library(gridExtra)

# enable parallel computation and allow STAN to automatically overwrite models
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE) 

# set working directory
## setwd("...")

## load data
data_EV <- read.csv("data/egg_viability_data.csv",sep=",")

# create a separate dataframe for each species
data_EV_Cmol <- data_EV %>% filter(species_id=="Cmol")
data_EV_Cpal <- data_EV %>% filter(species_id=="Cpal")
data_EV_Cqui <- data_EV %>% filter(species_id=="Cqui")
data_EV_Cthei <- data_EV %>% filter(species_id=="Cthei") 

## prepare input data for STAN model

# temp. points on which to calculate generated quantities
steps = 0.1
temp <- seq(0,45,steps) 
# a vector storing integer ids for species each observation belongs to
species_id = as.integer(as.factor(data_EV$species_id)) 
# a vector storing integer ids for experiment each observation belongs to
exp_id = as.integer(as.factor(data_EV$experiment_id)) 

# data list collecting all information needed to initialize the STAN model
data <- list(Nspecies = length(unique(data_EV$species_id)), # number of species
             Nexp = length(unique(data_EV$experiment_id)), # number of experiments
             Nobs = nrow(data_EV), # total number of observations
             x = data_EV$temperature, # data predictor (temperature)
             y = data_EV$trait, # data outcome (development rate)
             species_id = species_id, # species id of each observation
             exp_id = exp_id, # experiment id of each observation
             N_new = length(temp), # number of temp. points for generated quantities
             x_new = temp) # temp. points for generated quantities

## Fitting the STAN model or reading a fitted model

# fitting
#fit <- stan(file = 'stan_models/egg_viability_model.stan',
#            data = data, iter=4000, chains=4,
#.           control = list(adapt_delta=0.8, max_treedepth=12))

# save fitted model 
#saveRDS(fit, "model_fits/egg_viability_fit.rds")

# read fitted model
fit <- readRDS("model_fits/egg_viability_fit.rds")

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
f_new <- summary(fit, pars="f_new")$summary # model fit summary for each temperature (n=451) and each species (n=4) --> 1804 rows
f_new_mean <- matrix(f_new[,1], nrow=data$N_new) # mean of model fits
f_new_0025 <- matrix(f_new[,4], nrow=data$N_new) # 2.5% quantile of model fits
f_new_0975 <- matrix(f_new[,8], nrow=data$N_new) # 97.5% quantile of model fits

# Function to make basic plots of mean model fit (+ 95% credible interval) for each species
plot.species <- function(i, x, y, title){
  plot(temp,f_new_mean[,i],type="l",ylim=c(0,1),xlim=c(0,45),
       xlab="Temperature (°C)",
       ylab="Egg viability",
       lwd=2,cex.axis=1,cex.lab=1,main=title)
  lines(temp,f_new_0025[,i],lty=2,lwd=2)
  lines(temp,f_new_0975[,i],lty=2,lwd=2)
  points(x,y,cex=1,pch=19)
} 

# plot estimates of each species' expected response
#Cmol plot
plot.species(1, data_EV_Cmol$temperature,data_EV_Cmol$trait, title="Cules pip. molestus")

#Cpal plot
plot.species(2, data_EV_Cpal$temperature,data_EV_Cpal$trait, title="Cules pip. pallens")

#Cqui plot
plot.species(3, data_EV_Cqui$temperature,data_EV_Cqui$trait, title="Cules pip.quinque")

#Cthei plot
plot.species(4, data_EV_Cthei$temperature,data_EV_Cthei$trait, title="Cules theileri")

## Visualization of expected temperature response of a new species generated by sampling from hierachical prior
f_new_spec <- summary(fit, pars="f_new_spec")$summary # extract summary for expected response of a new species
plot(temp,f_new_spec[,1],type="l",ylim=c(0,1),ylab="Egg viability",
     lwd=2,cex.axis=1,cex.lab=1, main="new species")
lines(temp,f_new_spec[,4],lty=2,lwd=2)
lines(temp,f_new_spec[,8],lty=2,lwd=2)

### more sophisticated plots for article

## Plots of mean model fit (+ 95% credible interval) for each species

#Cmol plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,1], 
                     lowerCI = f_new_0025[,1], upperCI = f_new_0975[,1])

n_exp = length(unique(data_EV_Cmol$experiment_id))
n_total = nrow(data_EV_Cmol)

plot1 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "black", linewidth = 0.6) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "black", alpha = 0.15) +
  geom_point(data = data_EV_Cmol, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 0.9, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens molestus"))))+
  geom_richtext(aes(x = 32, y = 0.9), 
                label = paste("n<sub>exp</sub> =", n_exp,
                              "<br>n<sub>total</sub> =", n_total),
                size = 3, color = "black",
                fill = NA, label.color = NA,  # Transparent background
                hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),         
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 1))

#Cpal plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,2], 
                     lowerCI = f_new_0025[,2], upperCI = f_new_0975[,2])

n_exp = length(unique(data_EV_Cpal$experiment_id))
n_total = nrow(data_EV_Cpal)

plot2 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "black", linewidth = 0.6) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "black", alpha = 0.15) +
  geom_point(data = data_EV_Cpal, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 0.9, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens pallens"))))+
  geom_richtext(aes(x = 0, y = 0.9), 
                label = paste("n<sub>exp</sub> =", n_exp,
                              "<br>n<sub>total</sub> =", n_total),
                size = 3, color = "black",
                fill = NA, label.color = NA,  # Transparent background
                hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),         
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 1))

#Cqui plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,3], 
                     lowerCI = f_new_0025[,3], upperCI = f_new_0975[,3])

n_exp = length(unique(data_EV_Cqui$experiment_id))
n_total = nrow(data_EV_Cqui)

plot3 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "black", linewidth = 0.6) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "black", alpha = 0.15) +
  geom_point(data = data_EV_Cqui, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 0.9, shape = 19) +
  labs(title= expression(paste(italic("Culex quinquefasciatus"))))+
  geom_richtext(aes(x = 0, y = 0.9), 
                label = paste("n<sub>exp</sub> =", n_exp,
                              "<br>n<sub>total</sub> =", n_total),
                size = 3, color = "black",
                fill = NA, label.color = NA,  # Transparent background
                hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),         
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 1))

#Cthei plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,4], 
                     lowerCI = f_new_0025[,4], upperCI = f_new_0975[,4])

n_exp = length(unique(data_EV_Cthei$experiment_id))
n_total = nrow(data_EV_Cthei)

plot4 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "black", linewidth = 0.6) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "black", alpha = 0.15) +
  geom_point(data = data_EV_Cthei, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 0.9, shape = 19) +
  labs(title= expression(paste(italic("Culex theileri"))))+
  geom_richtext(aes(x = 0, y = 0.9), 
                label = paste("n<sub>exp</sub> =", n_exp,
                              "<br>n<sub>total</sub> =", n_total),
                size = 3, color = "black",
                fill = NA, label.color = NA,  # Transparent background
                hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),         
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 1))


# Plot three Culex species together (Figure 9)
plot_list = list(plot3, plot1, plot2)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2, label_size = 12,
                               align = "h", axis = "b", labels = c('A', 'B', 'C', 'D'))

y.grob <- textGrob(expression(paste("Egg viability")), 
                   gp=gpar(col="black", fontsize=10), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=10))

grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))

#ggsave("Figures/egg_viability.tiff", 
#       plot = grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob)),
#       width = 6, height = 4.6, 
#       dpi = 600, units = "in", compression = "lzw")
