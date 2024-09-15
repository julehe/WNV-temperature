### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to fit temperature-dependent functions for mosquito 
### susceptibility to WNV infection (infection probability) with a multi-study 
### Bayesian hierarchical model implemented in STAN 
###
### The script also includes the necessary code to reproduce Figures 11, SI4.3
###

## load libraries
library(rstan)
library(tidyverse)
library(mcmcplots)
library(bayesplot)
library(posterior)
library(grid)
library(gridExtra)
library(ggtext)

# enable parallel computation and allow STAN to automatically overwrite models
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE) 

# set working directory
## setwd("...")

## load data
data_infprob <- read.csv("data/infection_probability_data.csv") 

## prepare input data for STAN model

# temp. points on which to calculate generated quantities
steps = 0.1
temp <- seq(0,45,steps) 
# a vector storing integer ids for experiment each observation belongs to                           
exp_id = as.integer(as.factor(data_infprob$experiment_id))

# data list collecting all information needed to initialize the STAN model
data <- list(Nobs = nrow(data_infprob), # total number of observations
             Nexp = length(unique(data_infprob$experiment_id)), # number of experiments
             ns = as.integer(data_infprob$sample.size), # sample size for each observation
             x = data_infprob$temperature, # predictor (temperature)
             y = as.integer(data_infprob$nbr_positive), # outcome (percentage infected)
             exp_id = exp_id, # experiment id of each observation
             N_new = length(temp), # number of temp. points for generated quantities
             x_new = temp) # temp. points for generated quantities

# function that returns initial values for alpha so that chains do not start in 
# unplausible region
init_fun <- function(...) {
  list(alpha=runif(1,-10,-1)) 
}

## Fitting the STAN model or reading a fitted model

# fitting
#fit <- stan(file = 'stan_models/infection_probability_model.stan',
#           data = data, iter=4000,chains=4,
#           control = list(adapt_delta=0.99,max_treedepth=12), init=init_fun)

# save fitted model 
#saveRDS(fit, "model_fits/infection_probability_fit.rds")

# read fitted model
fit <- readRDS("model_fits/infection_probability_fit.rds")

## inspect the STAN model fit 

# print model summary
print(fit)

# print study level parameter estimates 
print(fit,pars=c("alpha"))
print(fit,pars=c("beta"))

# look at pairs plots to detect potential sampling problems
pairs(fit, pars=c("mu_alpha","mu_beta")) # population-level mean 
mcmc_pairs(fit, pars=c("mu_alpha","mu_beta"), 
           diag_fun="dens", off_diag_fun = "hex") # different visualization using the bayesplot package

# look at density- and traceplots to detect potential sampling problems
s = as.array(fit)
mcmc <- do.call(mcmc.list, plyr:::alply(s[, , -(length(s[1, 1, ]))], 2, as.mcmc))
denplot(mcmc, parms = c("mu_alpha","mu_beta"))
traplot(mcmc, parms = c("mu_alpha","mu_beta"))

### Visualizations

## population-level plot (expected temperature response in new experiment)
f_new <- summary(fit, pars="f_new_spec")$summary
f_new_mean <- matrix(f_new[,1],nrow=data$N_new) # mean of model fit
f_new_0025 <- matrix(f_new[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_0975 <- matrix(f_new[,8],nrow=data$N_new) # 97.5% quantile of model fit

temp_df = data.frame(temp = temp, mean = f_new_mean, 
                     lowerCI = f_new_0025, upperCI = f_new_0975)
plot_pop_level <- ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  labs(x= "", y= "",
       title = expression(paste("WNV in ",italic("Culex")," (population-level)"))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

plot_pop_level

## Plots of mean model fit (+ 95% credible interval) for each experiment

# NY99 in Cx. pipiens (Dohm et al.) 
f_new_Dohm <- summary(fit, pars="f_new_Dohm")$summary
f_new_Dohm_mean <- matrix(f_new_Dohm[,1],nrow=data$N_new) # mean of model fit
f_new_Dohm_0025 <- matrix(f_new_Dohm[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_Dohm_0975 <- matrix(f_new_Dohm[,8],nrow=data$N_new) # 97.5% quantile of model fit
temp_df = data.frame(temp = temp, mean = f_new_Dohm_mean, 
                     lowerCI = f_new_Dohm_0025, upperCI = f_new_Dohm_0975)

# extract data corresponding to this experiment
id = 1
data_temp <- data_infprob %>% filter(experiment_id %in% id)

# plot
plot1 <- ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_temp, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(x="",y="",
       title = expression(paste("WNV NY99 in ",italic("Culex pipiens"), " (Dohm)"))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

# WN02 in Cx. pipiens (Kilpatrick et al.) 
f_new_KPWN02 <- summary(fit, pars="f_new_KPWN02")$summary 
f_new_KPWN02_mean <- matrix(f_new_KPWN02[,1],nrow=data$N_new) # mean of model fit
f_new_KPWN02_0025 <- matrix(f_new_KPWN02[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_KPWN02_0975 <- matrix(f_new_KPWN02[,8],nrow=data$N_new) # 97.5% quantile of model fit
temp_df = data.frame(temp = temp, mean = f_new_KPWN02_mean, 
                     lowerCI = f_new_KPWN02_0025, upperCI = f_new_KPWN02_0975)

# extract data corresponding to this experiment
id = 2
data_temp <- data_infprob %>% filter(experiment_id %in% id)

# plot
plot2 <- ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_temp, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(x="",y="",
       title = expression(paste("WNV WN02 in ",italic("Culex pipiens"), " (Kilpatrick)"))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

# NY99 in Cx. pipiens (Kilpatrick et al.) 
f_new_KPNY99 <- summary(fit, pars="f_new_KPNY99")$summary
f_new_KPNY99_mean <- matrix(f_new_KPNY99[,1],nrow=data$N_new) # mean of model fit
f_new_KPNY99_0025 <- matrix(f_new_KPNY99[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_KPNY99_0975 <- matrix(f_new_KPNY99[,8],nrow=data$N_new) # 97.5% quantile of model fit
temp_df = data.frame(temp = temp, mean = f_new_KPNY99_mean, 
                     lowerCI = f_new_KPNY99_0025, upperCI = f_new_KPNY99_0975)

# extract data corresponding to this experiment
id = 3
data_temp <- data_infprob %>% filter(experiment_id %in% id)

# plot
plot3 <- ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_temp, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(x="",y="",
       title = expression(paste("WNV NY99 in ",italic("Culex pipiens"), " (Kilpatrick)"))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

# Plot all three experiments together (Figure 11)
plot_list = list(plot1, plot2, plot3)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B', 'C', 'D'))

y.grob <- textGrob(expression(paste("Mosquito infection probability")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/infection_probability.pdf", width=8.27, height=6.18)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()

### Comparison of models with different levels of between-experiment variability

# Fit or load alternative STAN model with increased between-experiment standard deviation
#fit <- stan(file = 'stan_models/infection_probability_model_increased_experiment_variability.stan',
#            data = data, iter=4000,chains=4,
#            control = list(adapt_delta=0.99, max_treedepth=12), 
#            init=init_fun)
#saveRDS(fit, "model_fits/infection_probability_fit_increased_experiment_variability.rds")
fit <- readRDS("model_fits/infection_probability_fit_increased_experiment_variability.rds")

# Extract population-level fit (expected temperature response in new experiment)
f_new <- summary(fit, pars="f_new_spec")$summary
f_new_mean <- matrix(f_new[,1],nrow=data$N_new) # mean of model fit
f_new_0025 <- matrix(f_new[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_0975 <- matrix(f_new[,8],nrow=data$N_new) # 97.5% quantile of model fit

# population-level plot (expected temperature response in new experiment) with 
# increased between-experiment variability
temp_df = data.frame(temp = temp, mean = f_new_mean, 
                     lowerCI = f_new_0025, upperCI = f_new_0975)
plot_pop_level_increased_var <- ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  labs(x= "", y= "",
       title = expression(paste("WNV in ",italic("Culex")," (population-level)"))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

# Extract fit and data for WN02 in Cx. pipiens (Kilpatrick et al.)
f_new_KPWN02 <- summary(fit, pars="f_new_KPWN02")$summary
f_new_KPWN02_mean <- matrix(f_new_KPWN02[,1],nrow=data$N_new) # mean of model fit
f_new_KPWN02_0025 <- matrix(f_new_KPWN02[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_KPWN02_0975 <- matrix(f_new_KPWN02[,8],nrow=data$N_new) # 97.5% quantile of model fit
temp_df = data.frame(temp = temp, mean = f_new_KPWN02_mean, 
                     lowerCI = f_new_KPWN02_0025, upperCI = f_new_KPWN02_0975)
id = 2
data_temp <- data_infprob %>% filter(experiment_id %in% id)

# plot for WN02 in Cx. pipiens (Kilpatrick et al.) with increased between-experiment 
# between-experiment variability
plot2_increased_var <- ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_temp, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(x="",y="",
       title = expression(paste("WNV WN02 in ",italic("Culex pipiens"), " (Kilpatrick)"))) +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

# Fit or load alternative STAN model with reduced between-experiment standard deviation
#fit <- stan(file = 'stan_models/infection_probability_model_reduced_experiment_variability.stan',
#            data = data, iter=4000,chains=4,
#            control = list(adapt_delta=0.99, max_treedepth=12), 
#            init=init_fun)
#saveRDS(fit, "model_fits/infection_probability_fit_reduced_experiment_variability.rds")
fit <- readRDS("model_fits/infection_probability_fit_reduced_experiment_variability.rds")

# Extract population-level fit (expected temperature response in new experiment)
f_new <- summary(fit, pars="f_new_spec")$summary
f_new_mean <- matrix(f_new[,1],nrow=data$N_new) # mean of model fit
f_new_0025 <- matrix(f_new[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_0975 <- matrix(f_new[,8],nrow=data$N_new) # 97.5% quantile of model fit

# population-level plot (expected temperature response in new experiment) with 
# reduced between-experiment variability
temp_df = data.frame(temp = temp, mean = f_new_mean, 
                     lowerCI = f_new_0025, upperCI = f_new_0975)
plot_pop_level_reduced_var <- ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  labs(x= "", y= "",
       title = expression(paste("WNV in ",italic("Culex")," (population-level)"))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))


# Extract fit and data for WN02 in Cx. pipiens (Kilpatrick et al.)
f_new_KPWN02 <- summary(fit, pars="f_new_KPWN02")$summary
f_new_KPWN02_mean <- matrix(f_new_KPWN02[,1],nrow=data$N_new) # mean of model fit
f_new_KPWN02_0025 <- matrix(f_new_KPWN02[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_KPWN02_0975 <- matrix(f_new_KPWN02[,8],nrow=data$N_new) # 97.5% quantile of model fit
temp_df = data.frame(temp = temp, mean = f_new_KPWN02_mean, 
                     lowerCI = f_new_KPWN02_0025, upperCI = f_new_KPWN02_0975)

id = 2
data_temp <- data_infprob %>% filter(experiment_id %in% id)

# plot for WN02 in Cx. pipiens (Kilpatrick et al.) with reduced between-experiment 
# between-experiment variability
plot2_reduced_var <- ggplot() + 
  geom_line(data = temp_df, aes(x = temp, y = mean), color = "red", linewidth = 0.8) +
  geom_ribbon(data = temp_df, 
              aes(x = temp, ymin = lowerCI, ymax = upperCI), fill = "red", alpha = 0.15) +
  geom_point(data = data_temp, aes(x = temperature, y = trait), size = 1.5, shape = 19) +
  labs(x="",y="",
       title = expression(paste("WNV WN02 in ",italic("Culex pipiens"), " (Kilpatrick)"))) +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 1))

# Combine model comparison plots together (Figure SI4.3)
plot_list = list(plot2_reduced_var, plot_pop_level_reduced_var,
                 plot2, plot_pop_level,
                 plot2_increased_var, plot_pop_level_increased_var)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B', 
                                                                   'C', 'D',
                                                                   'E', 'F'))

y.grob <- textGrob(expression(paste("Mosquito infection probability")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/infection_probability_different_levels_experiment_variability.pdf", width=8.27, height=9.27)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()


