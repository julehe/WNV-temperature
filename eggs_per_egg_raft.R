### Julian Heidecke, Heidelberg University, IWR 
### julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
###
### This script serves to plot data on the eggs per egg raft of Culex species
###
### The script includes the necessary code to reproduce Figure SI1.1
###

## load libraries
library(rstan)
library(tidyverse)
library(ggplot2)
library(ggtext)

## load data
data_eggs_per_raft <- read.csv("data/eggs_per_raft_data.csv",sep=",")

# create a separate dataframe for each species
data_eggs_per_raft_Cmol <- data_eggs_per_raft %>% filter(species_id=="Cmol")
data_eggs_per_raft_Cpal <- data_eggs_per_raft %>% filter(species_id=="Cpal")
data_eggs_per_raft_Cpip <- data_eggs_per_raft %>% filter(species_id=="Cpip")
data_eggs_per_raft_Cqui <- data_eggs_per_raft %>% filter(species_id=="Cqui")

# plot data (Figure SI1.1)

plot1 <- ggplot() + 
  geom_point(data = data_eggs_per_raft_Cmol, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 1.1, shape = 19) +
  labs(title= expression(paste(italic("Cx. pipiens molestus"))))+
  geom_richtext(aes(x = 0, y = 270), 
                label = paste("n<sub>exp</sub> =", length(unique(data_eggs_per_raft_Cmol$experiment_id)),
                              "<br>n<sub>total</sub> =", nrow(data_eggs_per_raft_Cmol)),
                size = 3, color = "black",
                fill = NA, label.color = NA,  # Transparent background
                hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  coord_cartesian(xlim = c(0,45),
                  ylim = c(0, 300))

plot2 <- ggplot() + 
  geom_point(data = data_eggs_per_raft_Cpal, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 1.1, shape = 19) +
  labs(title= expression(paste(italic("Cx. pipiens pallens"))))+
  geom_richtext(aes(x = 0, y = 270), 
                label = paste("n<sub>exp</sub> =", length(unique(data_eggs_per_raft_Cpal$experiment_id)),
                              "<br>n<sub>total</sub> =", nrow(data_eggs_per_raft_Cpal)),
                size = 3, color = "black",
                fill = NA, label.color = NA,  # Transparent background
                hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  coord_cartesian(xlim = c(0,45),
                  ylim = c(0, 300))

plot3 <- ggplot() + 
  geom_point(data = data_eggs_per_raft_Cpip, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 1.1, shape = 19) +
  labs(title= expression(paste(italic("Cx. pipiens"))))+
  geom_richtext(aes(x = 0, y = 270), 
                label = paste("n<sub>exp</sub> =", length(unique(data_eggs_per_raft_Cpip$experiment_id)),
                              "<br>n<sub>total</sub> =", nrow(data_eggs_per_raft_Cpip)),
                size = 3, color = "black",
                fill = NA, label.color = NA,  # Transparent background
                hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  coord_cartesian(xlim = c(0,45),
                  ylim = c(0, 300))

plot4 <- ggplot() + 
  geom_point(data = data_eggs_per_raft_Cqui, aes(x = temperature, y = trait, color=as.factor(experiment_id)), size = 1.1, shape = 19) +
  labs(title= expression(paste(italic("Cx. quinquefasciatus"))))+
  geom_richtext(aes(x = 0, y = 270), 
                label = paste("n<sub>exp</sub> =", length(unique(data_eggs_per_raft_Cqui$experiment_id)),
                              "<br>n<sub>total</sub> =", nrow(data_eggs_per_raft_Cqui)),
                size = 3, color = "black",
                fill = NA, label.color = NA,  # Transparent background
                hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  coord_cartesian(xlim = c(0,45),
                  ylim = c(0, 300))

# Plot four Culex species together (Figure SI1.1)
plot_list = list(plot3, plot4, plot1, plot2)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2, label_size = 12,
                               align = "h", axis = "b", labels = c('A', 'B', 'C', 'D'))

y.grob <- textGrob(expression(paste("Eggs per egg raft")), 
                   gp=gpar(col="black", fontsize=10), rot=90)

x.grob <- textGrob("Temperature (°C)", 
                   gp=gpar(col="black", fontsize=10))

grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))

#ggsave("Figures/eggs_per_raft.tiff", 
#       plot = grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob)),
#       width = 6, height = 5, 
#       dpi = 600, units = "in", compression = "lzw")
