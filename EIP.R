### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to fit temperature-dependent functions for the 
### extrinsic incubation period for WNV infection with a multi-study 
### Bayesian hierarchical model implemented in STAN 
###
### The script also includes the necessary code to reproduce Figures 12, SI4.4
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
data_EIP <- read.csv("data/EIP_data.csv")

# filter out data from studies that only looked at disseminated infections
data_EIP <- data_EIP %>% 
  filter(quantity == "percentage transmitting of infected" | quantity == "percentage transmitting of exposed")
# transform data that is given as "percentage transmitting of exposed" into
# "percentage transmitting of infected" by dividing the predicted "percentage infected"
data_EIP$trait[data_EIP$quantity == "percentage transmitting of exposed"] = data_EIP$trait[data_EIP$quantity == "percentage transmitting of exposed"]/data_EIP$prediction.infected[data_EIP$quantity == "percentage transmitting of exposed"]

## prepare input data for STAN model

# temp. points on which to calculate generated quantities
steps = 0.1
temp <- seq(0,45,steps) 
# days post infection on which to calculate generated quantities
dpi_new = seq(0,120,1)
# a vector storing integer ids for each experiment and temperature setting combination   
index_id = as.integer(as.factor(data_EIP$index))
# a vector storing the temperature setting of each experiment and temperature setting combination   
x_index = numeric(length(unique(index_id)))
for(i in 1:length(unique(index_id))){
  x_index[i]=data_EIP[data_EIP$index==i,]$temperature[1]
}
# a vector storing the experiment id of each experiment and temperature setting combination    
exp_index = numeric(length(unique(index_id)))
for(i in 1:length(unique(index_id))){
  exp_index[i]=data_EIP[data_EIP$index==i,]$experiment_id[1]
}

# data list collecting all information needed to initialize the STAN model
data <- list(Nobs = nrow(data_EIP), # total number of observations
             Nexp = length(unique(data_EIP$experiment_id)), # number of experiments
             Nindex = length(unique(index_id)), # number of experiment and temperature setting combinations
             dpi = data_EIP$dpi, # predictor (days post infection)
             id = index_id, # id of experiment and temperature setting combination
             y = data_EIP$trait, # outcome (percentage transmitting)
             N_new = length(temp), # number of temp. points for generated quantities
             x_new = temp, # temp. points for generated quantities
             N_dpi_new = length(dpi_new), # number of days post infection for generated quantities
             dpi_new = dpi_new, # days post infection for generated quantities
             x_index = x_index, # temperature setting of each experiment and temperature setting combination   
             exp_index = exp_index) # experiment id of each experiment and temperature setting combination   

## Fitting the STAN model or reading a fitted model

# fitting
#fit <- stan(file = 'stan_models/EIP_model.stan',
#           data = data, iter=4000,chains=4,
#           control = list(adapt_delta=0.99,max_treedepth=12))

# save fitted model 
#saveRDS(fit, "model_fits/EIP_fit.rds")

# read fitted model
fit <- readRDS("model_fits/EIP_fit.rds")

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

## population-level plot EIP vs temperature (expected temperature response in new experiment)
f_new_spec <- summary(fit, pars="f_new_spec")$summary
f_new_spec_mean <- matrix(f_new_spec[,1],nrow=data$N_new) # mean of model fit
f_new_spec_0025 <- matrix(f_new_spec[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_spec_0975 <- matrix(f_new_spec[,8],nrow=data$N_new) # 97.5% quantile of model fit

df_temp <- data.frame(tmp = temp, trait = f_new_spec_mean, 
                      lowerCI = f_new_spec_0025, upperCI = f_new_spec_0975)

plot_EIP_pop <- ggplot(data = df_temp, mapping = aes(x = tmp, y = trait,
                                                     ymin = lowerCI, 
                                                     ymax = upperCI)) +
  geom_line(linewidth = 0.8,color="red") +
  geom_ribbon(fill="red", alpha=0.15)+
  labs(title = "Population-level") + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 250)) 

plot_EIP_pop

## Plots of mean model fit (+ 95% credible interval) for each experiment

# Plots showing DPI vs transmitting mosquitoes 
# extract dpi predictions from model fit
p_new <- summary(fit, pars="p_new")$summary

# WN02 in Cx. pipiens (Kilpatrick et al.):

# extract data corresponding to this experiment
id = c(1,3,5,7)
data_temp <- data_EIP %>% filter(index %in% id)
# dpi predictions from model fit and bring into long format 
p_new_mean <- matrix(p_new[,1],nrow=data$N_dpi_new)[,id] # mean of model fit
p_new_mean <- matrix(p_new_mean, dimnames=list(t(outer(colnames(p_new_mean), rownames(p_new_mean), FUN=paste)), NULL))
p_new_0025 <- matrix(p_new[,4],nrow=data$N_dpi_new)[,id] # 2.5% quantile of model fit
p_new_0025 <- matrix(p_new_0025, dimnames=list(t(outer(colnames(p_new_0025), rownames(p_new_0025), FUN=paste)), NULL))
p_new_0975 <- matrix(p_new[,8],nrow=data$N_dpi_new)[,id] # 97.5% quantile of model fit
p_new_0975 <- matrix(p_new_0975, dimnames=list(t(outer(colnames(p_new_0975), rownames(p_new_0975), FUN=paste)), NULL))

# summarize data and predictions in dataframe and add temperature setting information
df_temp <- data.frame(dpi = dpi_new, trait = p_new_mean, tmp = as.factor(c(rep(c(32,22,18,15), each=data$N_dpi_new))), 
                      lowerCI = p_new_0025, upperCI = p_new_0975)
df_data_temp <- data.frame(x = data_temp$dpi, y = data_temp$trait, tmp = as.factor(data_temp$temperature))

# plot
plot1 <- ggplot(df_temp) +
  geom_line(mapping = aes(x = dpi, y = trait, group = tmp, color = tmp), linewidth = 0.8) +
  geom_ribbon(mapping = aes(x = dpi, ymin = lowerCI, ymax = upperCI, group = tmp, fill = tmp), alpha = 0.15,show.legend = FALSE) +
  geom_point(df_data_temp, mapping = aes(x=x, y=y, group = tmp, color = tmp)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = expression(paste("WNV WN02 in ",italic("Cx. pipiens")))) +  
  scale_color_discrete(labels = c("15°C", "18°C", "22°C", "32°C")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(.85,.725),
        legend.title = element_blank(),  
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), "pt"))


# NY99 in Cx. pipiens (Kilpatrick et al.):

# extract data corresponding to this experiment
id = c(2,4,6,8)
data_temp <- data_EIP %>% filter(index %in% id)
# dpi predictions from model fit and bring into long format 
p_new_mean <- matrix(p_new[,1],nrow=data$N_dpi_new)[,id] # mean of model fit
p_new_mean <- matrix(p_new_mean, dimnames=list(t(outer(colnames(p_new_mean), rownames(p_new_mean), FUN=paste)), NULL))
p_new_0025 <- matrix(p_new[,4],nrow=data$N_dpi_new)[,id] # 2.5% quantile of model fit
p_new_0025 <- matrix(p_new_0025, dimnames=list(t(outer(colnames(p_new_0025), rownames(p_new_0025), FUN=paste)), NULL))
p_new_0975 <- matrix(p_new[,8],nrow=data$N_dpi_new)[,id] # 97.5% quantile of model fit
p_new_0975 <- matrix(p_new_0975, dimnames=list(t(outer(colnames(p_new_0975), rownames(p_new_0975), FUN=paste)), NULL))

# summarize data and predictions in dataframe and add temperature setting information
df_temp <- data.frame(dpi = dpi_new, trait = p_new_mean, tmp = as.factor(c(rep(c(32,22,18,15), each=data$N_dpi_new))), 
                      lowerCI = p_new_0025, upperCI = p_new_0975)
df_data_temp <- data.frame(x = data_temp$dpi, y = data_temp$trait, tmp = as.factor(data_temp$temperature))

# plot
plot2 <- ggplot(df_temp) +
  geom_line(mapping = aes(x = dpi, y = trait, group = tmp, color = tmp), linewidth = 0.8) +
  geom_ribbon(mapping = aes(x = dpi, ymin = lowerCI, ymax = upperCI, group = tmp, fill = tmp), alpha = 0.15,show.legend = FALSE) +
  geom_point(df_data_temp, mapping = aes(x=x, y=y, group = tmp, color = tmp)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = expression(paste("WNV NY99 in ",italic("Cx. pipiens")))) +  
  scale_color_discrete(labels = c("15°C", "18°C", "22°C", "32°C")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position=c(.85,.725),
        legend.title = element_blank(),  
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), "pt")) 

# H442 in Cx. univittatus (Cornel et al.):

# extract data corresponding to this experiment
id = 9:12
data_temp <- data_EIP %>% filter(index %in% id)
# dpi predictions from model fit and bring into long format 
p_new_mean <- matrix(p_new[,1],nrow=data$N_dpi_new)[,id] # mean of model fit
p_new_mean <- matrix(p_new_mean, dimnames=list(t(outer(colnames(p_new_mean), rownames(p_new_mean), FUN=paste)), NULL))
p_new_0025 <- matrix(p_new[,4],nrow=data$N_dpi_new)[,id] # 2.5% quantile of model fit
p_new_0025 <- matrix(p_new_0025, dimnames=list(t(outer(colnames(p_new_0025), rownames(p_new_0025), FUN=paste)), NULL))
p_new_0975 <- matrix(p_new[,8],nrow=data$N_dpi_new)[,id] # 97.5% quantile of model fit
p_new_0975 <- matrix(p_new_0975, dimnames=list(t(outer(colnames(p_new_0975), rownames(p_new_0975), FUN=paste)), NULL))

# summarize data and predictions in dataframe and add temperature setting information
df_temp <- data.frame(dpi = dpi_new, trait = p_new_mean, tmp = as.factor(c(rep(c(30,26,18,14), each=data$N_dpi_new))), 
                      lowerCI = p_new_0025, upperCI = p_new_0975)
df_data_temp <- data.frame(x = data_temp$dpi, y = data_temp$trait, tmp = as.factor(data_temp$temperature))

# plot
plot3 <- ggplot(df_temp) +
  geom_line(mapping = aes(x = dpi, y = trait, group = tmp, color = tmp), linewidth = 0.8) +
  geom_ribbon(mapping = aes(x = dpi, ymin = lowerCI, ymax = upperCI, group = tmp, fill = tmp), alpha = 0.15,show.legend = FALSE) +
  geom_point(df_data_temp, mapping = aes(x=x, y=y, group = tmp, color=tmp)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = expression(paste("WNV H442 in ",italic("Cx. univittatus")))) +  
  scale_color_discrete(labels = c("14°C", "18°C", "26°C", "30°C")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
        #axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position=c(.85,.725),
        legend.title = element_blank(),  
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), "pt")) 

# NY99 in Cx. tarsalis (Reisen et al.):

# extract data corresponding to this experiment
id = 13:17
data_temp <- data_EIP %>% filter(index %in% id)
# dpi predictions from model fit and bring into long format 
p_new_mean <- matrix(p_new[,1],nrow=data$N_dpi_new)[,id] # mean of model fit
p_new_mean <- matrix(p_new_mean, dimnames=list(t(outer(colnames(p_new_mean), rownames(p_new_mean), FUN=paste)), NULL))
p_new_0025 <- matrix(p_new[,4],nrow=data$N_dpi_new)[,id] # 2.5% quantile of model fit
p_new_0025 <- matrix(p_new_0025, dimnames=list(t(outer(colnames(p_new_0025), rownames(p_new_0025), FUN=paste)), NULL))
p_new_0975 <- matrix(p_new[,8],nrow=data$N_dpi_new)[,id] # 97.5% quantile of model fit
p_new_0975 <- matrix(p_new_0975, dimnames=list(t(outer(colnames(p_new_0975), rownames(p_new_0975), FUN=paste)), NULL))

# summarize data and predictions in dataframe and add temperature setting information
df_temp <- data.frame(dpi = dpi_new, trait = p_new_mean, tmp = as.factor(c(rep(c(30,26,22,18,14), each=data$N_dpi_new))), 
                      lowerCI = p_new_0025, upperCI = p_new_0975)
df_data_temp <- data.frame(x = data_temp$dpi, y = data_temp$trait, tmp = as.factor(data_temp$temperature))

# plot
plot4 <- ggplot(df_temp) +
  geom_line(mapping = aes(x = dpi, y = trait, group = tmp, color = tmp), linewidth = 0.8) +
  geom_ribbon(mapping = aes(x = dpi, ymin = lowerCI, ymax = upperCI, group = tmp, fill = tmp), alpha = 0.15,show.legend = FALSE) +
  geom_point(df_data_temp, mapping = aes(x=x, y=y, group = tmp, color = tmp)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = expression(paste("WNV NY99 in ",italic("Cx. tarsalis")))) +  
  scale_color_discrete(labels = c("14°C", "18°C", "22°C", "26°C", "30°C")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
        #axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position=c(.85,.25),
        legend.title = element_blank(),  
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), "pt")) 

# Combine plots showing DPI vs transmitting mosquitoes (part of Figure 12)
plot_list = list(plot1, plot2, plot3, plot4)
plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2, 
                               align = "h", axis = "b", labels = c('A','B','C','D'))

y.grob <- textGrob(expression(paste("Percentage of transmitting mosquitoes")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Days post infection (dpi)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/EIP1.pdf", width=8.3, height=7)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()

# Plots showing EIP vs temperature at the experiment-level
f_new <- summary(fit, pars="f_new")$summary 
f_new_mean <- matrix(f_new[,1],nrow=data$N_new) # mean of model fit
f_new_0025 <- matrix(f_new[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_0975 <- matrix(f_new[,8],nrow=data$N_new) # 97.5% quantile of model fit

# prepare dataframe storing EIP predictions for each experiment and add experiment identifier column
df_temp <- data.frame(tmp = temp, trait = c(f_new_mean[,1],f_new_mean[,2],f_new_mean[,3],f_new_mean[,4]),
                      lowerCI = c(f_new_0025[,1],f_new_0025[,2],f_new_0025[,3],f_new_0025[,4]),
                      upperCI = c(f_new_0975[,1],f_new_0975[,2],f_new_0975[,3],f_new_0975[,4]),
                      experiment = c(rep("1pip_WN02",length(temp)),rep("2pip_NY99",length(temp)),
                                     rep("3uni_H442",length(temp)),rep("4tar_NY99",length(temp))))

# plot EIP vs temperature for WNV02 in Cx. pipiens (Kilpatrick et al.)
plot_EIP_WN02Cpip <- ggplot(data = filter(df_temp,experiment == "1pip_WN02"), 
                                mapping = aes(x = tmp, y = trait,
                                              ymin = lowerCI, ymax = upperCI)) +
  geom_line(linewidth = 0.8,color="red") +
  geom_ribbon(fill="red", alpha=0.15)+
  labs(y = "", x = "", title = expression(paste("WNV WN02 in ",italic("Cx. pipiens")))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position=c(.7,.8),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 250)) 

# plot EIP vs temperature for all experiments
plot_EIP_exp <- ggplot(data = df_temp, mapping = aes(x = tmp, y = trait,
                                                     color = experiment,
                                                     group = experiment, 
                                                     ymin = lowerCI, 
                                                     ymax = upperCI,
                                                     fill = experiment)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(color=NA, alpha=0.15)+
  scale_color_manual(values = c("black","orange","blue","green4"),
                     labels = c(expression(paste("WNV WN02 in ",italic("Cx. pipiens"))),
                                expression(paste("WNV NY99 in ",italic("Cx. pipiens"))),
                                expression(paste("WNV H442 in ",italic("Cx. univittatus"))),
                                expression(paste("WNV NY99 in ",italic("Cx. tarsalis"))))) +
  scale_fill_manual(values = c("black","orange","blue","green4"),
                    guide = "none") +
  labs(title = "", x = "Temperature (°C)") + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.position=c(.7,.8),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 250)) 

# prepare plot EIP vs temperature for all experiments as plot grid to fit format of DPI plots (part of Figure 12)
plot_list = list(plot_EIP_exp)
plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2, 
                               align = "h", axis = "b", labels = c('E'))

y.grob <- textGrob(expression(paste("Extrinsic incubation period (days)")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

#pdf("Figures/EIP2.pdf", width=8.3, height=4)
grid.arrange(arrangeGrob(plot_grid, left = y.grob))
#dev.off()

### Comparison of models with different levels of between-experiment variability

# Fit or load alternative STAN model with increased between-experiment standard deviation

#fit <- stan(file = 'stan_models/EIP_model_increased_experiment_variability.stan', 
#            data = data, iter=4000, chains=4, 
#            control = list(adapt_delta=0.8, max_treedepth = 12))

#saveRDS(fit, "model_fits/EIP_fit_increased_experiment_variability.rds")
fit <- readRDS("model_fits/EIP_fit_increased_experiment_variability.rds")

# Extract population-level fit (expected temperature response of EIP in new experiment)
f_new_spec <- summary(fit, pars="f_new_spec")$summary
f_new_spec_mean <- matrix(f_new_spec[,1],nrow=data$N_new) # mean of model fit
f_new_spec_0025 <- matrix(f_new_spec[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_spec_0975 <- matrix(f_new_spec[,8],nrow=data$N_new) # 97.5% quantile of model fit

# population-level plot (expected temperature response in new experiment) with 
# increased between-experiment variability
df_temp <- data.frame(tmp = temp, trait = f_new_spec_mean, 
                      lowerCI = f_new_spec_0025, upperCI = f_new_spec_0975)
plot_EIP_pop_increased_var <- ggplot(data = df_temp, mapping = aes(x = tmp, y = trait,
                                                     ymin = lowerCI, 
                                                     ymax = upperCI)) +
  geom_line(linewidth = 0.8,color="red") +
  geom_ribbon(fill="red", alpha=0.15)+
  labs(title = "Population-level") + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 250)) 

# Extract fit and data for WN02 in Cx. pipiens (Kilpatrick et al.) 
f_new <- summary(fit, pars="f_new")$summary 
f_new_mean <- matrix(f_new[,1],nrow=data$N_new) # mean of model fit
f_new_0025 <- matrix(f_new[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_0975 <- matrix(f_new[,8],nrow=data$N_new) # 97.5% quantile of model fit

# plot for WN02 in Cx. pipiens (Kilpatrick et al.) with increased between-experiment 
# between-experiment variability
df_temp <- data.frame(tmp = temp, trait = c(f_new_mean[,1],f_new_mean[,2],f_new_mean[,3],f_new_mean[,4]),
                      lowerCI = c(f_new_0025[,1],f_new_0025[,2],f_new_0025[,3],f_new_0025[,4]),
                      upperCI = c(f_new_0975[,1],f_new_0975[,2],f_new_0975[,3],f_new_0975[,4]),
                      experiment = c(rep("1pip_WN02",length(temp)),rep("2pip_NY99",length(temp)),
                                     rep("3uni_H442",length(temp)),rep("4tar_NY99",length(temp))))

plot_EIP_WN02Cpip_increased_var <- ggplot(data = filter(df_temp,experiment == "1pip_WN02"), 
                                mapping = aes(x = tmp, y = trait,
                                              ymin = lowerCI, ymax = upperCI)) +
  geom_line(linewidth = 0.8,color="red") +
  geom_ribbon(fill="red", alpha=0.15)+
  labs(y = "", x = "", title = expression(paste("WNV WN02 in ",italic("Cx. pipiens")))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position=c(.7,.8),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 250)) 

# Fit or load alternative STAN model with reduced between-experiment standard deviation

#fit <- stan(file = 'stan_models/EIP_model_reduced_experiment_variability.stan', 
#            data = data, iter=4000, chains=4, 
#            control = list(adapt_delta=0.8, max_treedepth = 12))

#saveRDS(fit, "model_fits/EIP_fit_reduced_experiment_variability.rds")
fit <- readRDS("model_fits/EIP_fit_reduced_experiment_variability.rds")

# Extract population-level fit (expected temperature response of EIP in new experiment)
f_new_spec <- summary(fit, pars="f_new_spec")$summary
f_new_spec_mean <- matrix(f_new_spec[,1],nrow=data$N_new) # mean of model fit
f_new_spec_0025 <- matrix(f_new_spec[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_spec_0975 <- matrix(f_new_spec[,8],nrow=data$N_new) # 97.5% quantile of model fit

# population-level plot (expected temperature response in new experiment) with 
# reduced between-experiment variability
df_temp <- data.frame(tmp = temp, trait = f_new_spec_mean, 
                      lowerCI = f_new_spec_0025, upperCI = f_new_spec_0975)
plot_EIP_pop_reduced_var <- ggplot(data = df_temp, mapping = aes(x = tmp, y = trait,
                                                                   ymin = lowerCI, 
                                                                   ymax = upperCI)) +
  geom_line(linewidth = 0.8,color="red") +
  geom_ribbon(fill="red", alpha=0.15)+
  labs(title = "Population-level") + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 250)) 

# Extract fit and data for WN02 in Cx. pipiens (Kilpatrick et al.) 
f_new <- summary(fit, pars="f_new")$summary 
f_new_mean <- matrix(f_new[,1],nrow=data$N_new) # mean of model fit
f_new_0025 <- matrix(f_new[,4],nrow=data$N_new) # 2.5% quantile of model fit
f_new_0975 <- matrix(f_new[,8],nrow=data$N_new) # 97.5% quantile of model fit

# plot for WN02 in Cx. pipiens (Kilpatrick et al.) with reduced between-experiment 
# between-experiment variability
df_temp <- data.frame(tmp = temp, trait = c(f_new_mean[,1],f_new_mean[,2],f_new_mean[,3],f_new_mean[,4]),
                      lowerCI = c(f_new_0025[,1],f_new_0025[,2],f_new_0025[,3],f_new_0025[,4]),
                      upperCI = c(f_new_0975[,1],f_new_0975[,2],f_new_0975[,3],f_new_0975[,4]),
                      experiment = c(rep("1pip_WN02",length(temp)),rep("2pip_NY99",length(temp)),
                                     rep("3uni_H442",length(temp)),rep("4tar_NY99",length(temp))))

plot_EIP_WN02Cpip_reduced_var <- ggplot(data = filter(df_temp,experiment == "1pip_WN02"), 
                                          mapping = aes(x = tmp, y = trait,
                                                        ymin = lowerCI, ymax = upperCI)) +
  geom_line(linewidth = 0.8,color="red") +
  geom_ribbon(fill="red", alpha=0.15)+
  labs(y = "", x = "", title = expression(paste("WNV WN02 in ",italic("Cx. pipiens")))) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position=c(.7,.8),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 250)) 

# Combine model comparison plots together (Figure SI4.4)
plot_list = list(plot_EIP_WN02Cpip_reduced_var, plot_EIP_pop_reduced_var,
                 plot_EIP_WN02Cpip, plot_EIP_pop,
                 plot_EIP_WN02Cpip_increased_var, plot_EIP_pop_increased_var)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2,
                               align = "h", axis = "b", labels = c('A', 'B', 
                                                                   'C', 'D',
                                                                   'E', 'F'))

y.grob <- textGrob(expression(paste("Extrinsic incubation period (days)")), 
                   gp=gpar(col="black", fontsize=14), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=14))

#pdf("Figures/EIP_different_levels_experiment_variability.pdf", width=8.27, height=9.27)
grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))
#dev.off()

