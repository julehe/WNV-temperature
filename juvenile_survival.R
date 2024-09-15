### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to fit temperature-dependent functions for mosquito juvenile
### survival (from egg hatch to adult emergence) using a multi-species 
### Bayesian hierarchical model implemented in STAN 
###
### The script also includes the necessary code to reproduce Figure 7
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
data_sur <- read.csv("data/juvenile_survival_data.csv",sep=",")

# create a separate dataframe for each species
data_sur_Anig <- data_sur %>% filter(species_id=="Anig")
data_sur_Asol <- data_sur %>% filter(species_id=="Asol")
data_sur_Atri <- data_sur %>% filter(species_id=="Atri")
data_sur_Avex <- data_sur %>% filter(species_id=="Avex")
data_sur_Cino <- data_sur %>% filter(species_id=="Cino")
data_sur_Cmel <- data_sur %>% filter(species_id=="Cmel")
data_sur_Cmol <- data_sur %>% filter(species_id=="Cmol")
data_sur_Cpal <- data_sur %>% filter(species_id=="Cpal")
data_sur_Cpip <- data_sur %>% filter(species_id=="Cpip")
data_sur_Cqui <- data_sur %>% filter(species_id=="Cqui")
data_sur_Cres <- data_sur %>% filter(species_id=="Cres") 
data_sur_Csal <- data_sur %>% filter(species_id=="Csal")
data_sur_Ctar <- data_sur %>% filter(species_id=="Ctar")

## prepare input data for STAN model

# temp. points on which to calculate generated quantities
steps = 0.1
temp <- seq(0,45,steps) 
# a vector storing integer ids for species each observation belongs to
species_id = as.integer(as.factor(data_sur$species_id)) 
# a vector storing integer ids for experiment each observation belongs to
exp_id = as.integer(as.factor(data_sur$experiment_id)) 

# data list collecting all information needed to initialize the STAN model
data <- list(Nspecies = length(unique(data_sur$species_id)), # number of species
             Nexp = length(unique(data_sur$experiment_id)), # number of experiments
             Nobs = nrow(data_sur), # total number of observations
             x = data_sur$temperature, # predictor (temperature)
             y = data_sur$trait, # outcome (survival rate)
             species_id = species_id, # species id of each observation
             exp_id = exp_id, # experiment id of each observation
             N_new = length(temp), # number of temp. points for generated quantities
             x_new = temp) # temp. points for generated quantities

## Fitting the STAN model or reading a fitted model

# fitting
#fit <- stan(file = 'stan_models/juvenile_survival_model.stan',
#            data = data, iter=4000, chains=4, 
#            control = list(adapt_delta=0.8, max_treedepth=12))

#fit <- stan(file = 'code/larva to adult survival (quadratic)/independent priors (ip)/study effect (se)/homosced error (hom)/surLA_ip_se_hom_test_qlogn.stan',
#            data = data, iter=4000, chains=4, control = list(adapt_delta=0.99, max_treedepth=12))

# save fitted model 
#saveRDS(fit, "model_fits/juvenile_survival_fit.rds")

# read fitted model
fit <- readRDS("model_fits/juvenile_survival_fit.rds")

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
f_new <- summary(fit, pars="f_new")$summary # model fit summary for each temperature (n=451) and each species (n=13) --> 5863 rows
f_new_mean <- matrix(f_new[,1],nrow=data$N_new) # mean of f_new
f_new_0025 <- matrix(f_new[,4],nrow=data$N_new) # 2.5% of f_new
f_new_0975 <- matrix(f_new[,8],nrow=data$N_new) # 97.5% of f_new

# Function to make basic plots of mean model fit (+ 95% credible interval) for each species
plot.species <- function(i, x, y, title){
  plot(temp,f_new_mean[,i],type="l",ylim=c(0,1),
       xlab="Temperature (°C)",
       ylab="Larva to Adult survival",
       lwd=2,cex.axis=1,cex.lab=1, main=title)
  lines(temp,f_new_0025[,i],lty=2,lwd=2)
  lines(temp,f_new_0975[,i],lty=2,lwd=2)
  points(x,y,cex=1,pch=19)
} 

# plot estimates of each species' expected response
#Anig plot
plot.species(1, data_sur_Anig$temperature,data_sur_Anig$trait, title="Aedes nigromaculis")

#Asol plot
plot.species(2, data_sur_Asol$temperature,data_sur_Asol$trait, title="Aedes solicitans")

#Atri plot
plot.species(3, data_sur_Atri$temperature,data_sur_Atri$trait, title="Aedes triseratus")

#Avex plot
plot.species(4, data_sur_Avex$temperature,data_sur_Avex$trait, title="Aedes vexans")

#Cino plot
plot.species(5, data_sur_Cino$temperature,data_sur_Cino$trait, title="Culiseta inornata")

#Cmel plot
plot.species(6, data_sur_Cmel$temperature,data_sur_Cmel$trait, title="Culiseta melanura")

#Cmol plot
plot.species(7, data_sur_Cmol$temperature,data_sur_Cmol$trait, title="Culex pip. molestus")

#Cpal plot
plot.species(8, data_sur_Cpal$temperature,data_sur_Cpal$trait, title="Culex pip. pallens")

#Cpip plot
plot.species(9, data_sur_Cpip$temperature,data_sur_Cpip$trait, title="Culex pipiens")

#Cqui plot
plot.species(10, data_sur_Cqui$temperature,data_sur_Cqui$trait, title="Culex quinque")

#Cres plot
plot.species(11, data_sur_Cres$temperature,data_sur_Cres$trait, title="Culex restuans")

#Csal plot
plot.species(12, data_sur_Csal$temperature,data_sur_Csal$trait, title="Culex salinarius")

#Ctar plot
plot.species(13, data_sur_Ctar$temperature,data_sur_Ctar$trait, title="Culex tarsalis")

## Visualization of expected temperature response of a new species generated by sampling from hierachical prior
f_new_spec <- summary(fit, pars="f_new_spec")$summary # extract summary for expected response of a new species
plot(temp,f_new_spec[,1],type="l",ylim=c(0,1),ylab="Larva to Adult survival",
     lwd=2,cex.axis=1,cex.lab=1, main="new species")
lines(temp,f_new_spec[,4],lty=2,lwd=2)
lines(temp,f_new_spec[,8],lty=2,lwd=2)

### more sophisticated plots for article

## Plots of mean model fit (+ 95% credible interval) for each species

#Anig plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,1], 
                     lowerCI = f_new_0025[,1], upperCI = f_new_0975[,1])
plot1 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Anig, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Aedes nigromaculis"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Asol plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,2], 
                     lowerCI = f_new_0025[,2], upperCI = f_new_0975[,2])
plot2 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Asol, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Aedes solicitans"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Atri plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,3], 
                     lowerCI = f_new_0025[,3], upperCI = f_new_0975[,3])
plot3 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Atri, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Aedes triseriatus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Avex plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,4], 
                     lowerCI = f_new_0025[,4], upperCI = f_new_0975[,4])
plot4 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Avex, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Aedes vexans"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Cino plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,5], 
                     lowerCI = f_new_0025[,5], upperCI = f_new_0975[,5])
plot5 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Cino, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culiseta inornata"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Cmel plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,6], 
                     lowerCI = f_new_0025[,6], upperCI = f_new_0975[,6])
plot6 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Cmel, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culiseta melanura"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Cmol plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,7], 
                     lowerCI = f_new_0025[,7], upperCI = f_new_0975[,7])
plot7 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Cmol, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens molestus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Cpal plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,8], 
                     lowerCI = f_new_0025[,8], upperCI = f_new_0975[,8])
plot8 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Cpal, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens pallens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Cpip plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,9], 
                     lowerCI = f_new_0025[,9], upperCI = f_new_0975[,9])
plot9 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Cpip, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex pipiens"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Cqui plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,10], 
                     lowerCI = f_new_0025[,10], upperCI = f_new_0975[,10])
plot10 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Cqui, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex quinquefasciatus"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Cres plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,11], 
                     lowerCI = f_new_0025[,11], upperCI = f_new_0975[,11])
plot11 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Cres, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex restuans"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

#Csal plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,12], 
                     lowerCI = f_new_0025[,12], upperCI = f_new_0975[,12])
plot12 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Csal, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex salinarius"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1)) 

#Ctar plot
temp_df = data.frame(temp = temp, mean = f_new_mean[,13], 
                     lowerCI = f_new_0025[,13], upperCI = f_new_0975[,13])
plot13 = ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_sur_Ctar, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(title= expression(paste(italic("Culex tarsalis"))))+
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1)) 

# Plot six Culex species together (Figure 7)
plot_list = list(plot9, plot10, plot7, plot8, plot11, plot13)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B', 'C', 'D', 'E', 'F'))

y.grob <- textGrob(expression(paste("Larva to adult survival")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/juvenile_survival.pdf", width=8.27, height=9.27)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()

