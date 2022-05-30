##### corHMM Plots #####
# load packages
library(ggplot2)
library(cowplot)

# load data
data <- read.csv(file="corhmm_group.csv", header=TRUE)

# load colours
colour_1 <- c("#e41a1c",
              "#377eb8",
              "#4daf4a",
              "#984ea3",
              "#ff7f00",
              "#f781bf","#a65628","#feb24c")

colour_2 <- c("#1b9e77","#e6ab02",
              "#d95f02",
              "#7570b3",
              "#e7298a",
              "#66a61e",
              "#a6761d",
              "#666666")

pa <- ggplot() + geom_point(aes(x=seq_along(data$global), y=data$global, colour= data$Mutualism, shape=as.factor(data$sig), size=0.0001)) +
  ylab("Weighted Difference in Diversification Rate") + xlab("Datasets") + theme_classic() + scale_color_manual(values = colour_2) + scale_shape_manual(values = c(21,19)) +ylim(-1,2)
pa <- pa + geom_hline(yintercept=0, linetype="dotted") +  geom_point(size = 5)

pb <- ggplot() + geom_point(aes(x=seq_along(data$global), y=data$global, colour= data$Class, shape=as.factor(data$sig), size=0.0001)) +
  ylab("Weighted Difference in Diversification Rate") + xlab("Datasets") + theme_classic() +  scale_color_manual(values = colour_1) + scale_shape_manual(values = c(21,19)) +ylim(-1,2)
pb <- pb + geom_hline(yintercept=0, linetype="dotted") +  geom_point(size = 5)




