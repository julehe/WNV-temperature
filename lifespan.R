### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to fit temperature-dependent functions for mosquito 
### adult lifespan using a multi-species Bayesian hierarchical model 
### implemented in STAN 
###
### The script also includes the necessary code to reproduce Figure 8
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
data_lf <- read.csv("data/lifespan_data.csv",sep=",")

# create a separate dataframe for each species
data_lf_Cpip <- data_lf %>% filter(species_id=="Cpip")
data_lf_Cqui <- data_lf %>% filter(species_id=="Cqui")
data_lf_Ctar <- data_lf %>% filter(species_id=="Ctar")
data_lf_Cres <- data_lf %>% filter(species_id=="Cres")
data_lf_Cpal <- data_lf %>% filter(species_id=="Cpal")
data_lf_Atae <- data_lf %>% filter(species_id=="Atae")
data_lf_Cmol <- data_lf %>% filter(species_id=="Cmol")

## prepare input data for STAN model

# temp. points on which to calculate generated quantities
steps = 0.1
temp <- seq(0,45,steps) # temp. points on which to calculate generated quantities
# a vector storing integer ids for species each observation belongs to
species_id = as.integer(as.factor(data_lf$species_id)) 
# a vector storing integer ids for experiment each observation belongs to
exp_id = as.integer(as.factor(data_lf$experiment_id)) 

# data list collecting all information needed to initialize the STAN model
data <- list(Nspecies = length(unique(data_lf$species_id)), # number of species
             Nexp = length(unique(data_lf$experiment_id)), # number of experiments
             Nobs = nrow(data_lf), # total number of observations
             x = data_lf$temperature, # predictor (temperature)
             y = data_lf$trait, # outcome (lifespan)
             species_id = species_id, # species id of each observation
             exp_id = exp_id, # experiment id of each observation
             N_new = length(temp), # number of temp. points for generated quantities
             x_new = temp) # temp. points for generated quantities

# function that returns initial values for Tmax so that chains do not start in 
# unplausible region
init_fun <- function(...) {
  list(mu_Tmax=runif(1,30,40))
}

## Fitting the STAN model or reading a fitted model

# fitting
#fit <- stan(file = 'stan_models/lifespan_model.stan',
#            data = data, iter=4000,chains=4,
#.           control = list(adapt_delta=0.99, max_treedepth=12),init=init_fun)

# save fitted model 
#saveRDS(fit, "model_fits/lifespan_fit.rds")

# read fitted model
fit <- readRDS("model_fits/lifespan_fit.rds")

## inspect STAN model fit 
# model summary
print(fit)

# print between species and experiment variability estimates 
print(fit,pars=c("sigma_exp_Tmax","sigma_exp_beta"),probs = c(0.025,0.975))
print(fit,pars=c("sigma_Tmax","sigma_beta"),probs = c(0.025,0.975))

# print species level parameter estimates 
print(fit,pars=c("beta"),probs = c(0.025,0.975))
print(fit,pars=c("Tmax"),probs = c(0.025,0.975))

# look at pairs plots to detect potential sampling problems
pairs(fit, pars=c("mu_beta", "mu_Tmax",
                  "sigma_beta","sigma_Tmax",
                  "sigma_exp_beta","sigma_exp_Tmax")) # population-level mean and between species/experiment varibility

mcmc_pairs(fit, pars=c("sigma_exp_beta","sigma_exp_Tmax"), 
           diag_fun="dens", off_diag_fun = "hex") # different visualization using the bayesplot package

# look at density- and traceplots to detect potential sampling problems
s = as.array(fit)
mcmc <- do.call(mcmc.list, plyr:::alply(s[, , -(length(s[1, 1, ]))], 2, as.mcmc))
denplot(mcmc, parms = c("mu_Tmax","mu_beta","sigma_Tmax","sigma_beta"))
denplot(mcmc, parms = c("sigma_exp")) 
traplot(mcmc, parms = c("mu_beta","sigma_beta","mu_Tmax","sigma_Tmax")) 
traplot(mcmc, parms = c("sigma_exp")) 

### Visualizations

## Basic plots of mean model fit (+ 95% credible interval) for each species

# Extract model fit summary (generated quantity)
f_new <- summary(fit, pars="f_new")$summary # model fit summary for each temperature (n=451) and each species (n=7) --> 3157 rows
f_new_mean <- matrix(f_new[,1],nrow=data$N_new) # mean of f_new
f_new_0025 <- matrix(f_new[,4],nrow=data$N_new) # 2.5% of f_new
f_new_0975 <- matrix(f_new[,8],nrow=data$N_new) # 97.5% of f_new

# introduce cut off at lowest observed temperature (14°C)
for(i in 1:(14/steps)){
  f_new_mean[i,] <- f_new_mean[(14/steps+1),]
  f_new_0025[i,] <- f_new_0025[(14/steps+1),]
  f_new_0975[i,] <- f_new_0975[(14/steps+1),]
}


# Function to make basic plots of mean model fit (+ 95% credible interval) for each species
plot.species <- function(i, x, y, title){
  plot(temp,f_new_mean[,i],type="l",ylim=c(0,160),
       xlab="Temperature (°C)",
       ylab="Adult lifespan",
       lwd=2,cex.axis=1,cex.lab=1, main=title)
  lines(temp,f_new_0025[,i],lty=2,lwd=2)
  lines(temp,f_new_0975[,i],lty=2,lwd=2)
  points(x,y,cex=1,pch=19)
} 

# plot estimates of each species' expected response
#Atae plot
plot.species(1, data_lf_Atae$temperature,data_lf_Atae$trait, title="Aedes taeniorhynchus")

#Cmol plot
plot.species(2, data_lf_Cmol$temperature,data_lf_Cmol$trait, title="Culex pip. molestus")

#Cpal plot
plot.species(3, data_lf_Cpal$temperature,data_lf_Cpal$trait, title="Culex pip. pallens")

#Cpip plot
plot.species(4, data_lf_Cpip$temperature,data_lf_Cpip$trait, title="Culex pipiens")

#Cqui plot
plot.species(5, data_lf_Cqui$temperature,data_lf_Cqui$trait, title="Culex quinque")

#Cres plot
plot.species(6, data_lf_Cres$temperature,data_lf_Cres$trait, title="Culex restuans")

#Ctar plot
plot.species(7, data_lf_Ctar$temperature,data_lf_Ctar$trait, title="Culex tarsalis")

## Visualization of expected temperature response of a new species generated by sampling from hierachical prior
f_new_spec <- summary(fit, pars="f_new_spec")$summary # extract summary for expected response of a new species
f_new_spec_mean <- matrix(f_new_spec[,1],nrow=data$N_new) # mean of f_new_spec
f_new_spec_0025 <- matrix(f_new_spec[,4],nrow=data$N_new) # 2.5% of f_new_spec
f_new_spec_0975 <- matrix(f_new_spec[,8],nrow=data$N_new) # 97.5% of f_new_spec

# introduce cut off at lowest observed temperature (14°C)
for(i in 1:(14/steps)){
  f_new_spec_mean[i,] <- f_new_spec_mean[(14/steps+1),]
  f_new_spec_0025[i,] <- f_new_spec_0025[(14/steps+1),]
  f_new_spec_0975[i,] <- f_new_spec_0975[(14/steps+1),]
}

plot(temp,f_new_spec_mean,type="l",ylim=c(0,250),ylab="Adult lifespan",
     lwd=2,cex.axis=1,cex.lab=1, main="new species")
lines(temp,f_new_spec_0025,lty=2,lwd=2)
lines(temp,f_new_spec_0975,lty=2,lwd=2)

### more sophisticated plots for article

## Plots of mean model fit (+ 95% credible interval) for each species

#Atae plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,1], 
                     lowerCI = f_new_0025[,1], upperCI = f_new_0975[,1])
plot1 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_lf_Atae, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Aedes taeniorhynchus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 160))

#Cmol plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,2], 
                     lowerCI = f_new_0025[,2], upperCI = f_new_0975[,2])
plot2 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_lf_Cmol, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens molestus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 160))

#Cpal plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,3], 
                     lowerCI = f_new_0025[,3], upperCI = f_new_0975[,3])
plot3 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_lf_Cpal, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens pallens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 160))

#Cpip plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,4], 
                     lowerCI = f_new_0025[,4], upperCI = f_new_0975[,4])
plot4 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_lf_Cpip, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 160))

#Cqui plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,5], 
                     lowerCI = f_new_0025[,5], upperCI = f_new_0975[,5])
plot5 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_lf_Cqui, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex quinquefasciatus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 160))

#Cres plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,6], 
                     lowerCI = f_new_0025[,6], upperCI = f_new_0975[,6])
plot6 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_lf_Cres, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex restuans"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 160))

#Ctar plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,7], 
                     lowerCI = f_new_0025[,7], upperCI = f_new_0975[,7])
plot7 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_lf_Ctar, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex tarsalis"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 160)) 

# Plot six Culex species together (Figure 8)
plot_list = list(plot4, plot5, plot2, plot3, plot6, plot7)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B', 'C', 'D', 'E', 'F'))

y.grob <- textGrob(expression(paste("Adult lifespan (days)")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

pdf("Figures/lifespan.pdf", width=8.27, height=9.27)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
dev.off()


