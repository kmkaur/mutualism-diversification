##### RevBayes Weighted Diversification Rates #####

#l oad packages
library(ggplot2)
library(dplyr)
library(cowplot)
library(phytools)

# load RevBayes log files with speciation and extinction rates
files <- list.files(pattern="*.log")
data <- lapply(files, read.table, header=TRUE)

# make empty lists
dat_div_split <- list()
start <- list()
end <- list()
HiSSE_types <- list()
dat_div <- list()
dat_div_split <- list()
A1_0 <- list()
B1_0 <- list()
A1_0df <- list()
B1_0df <- list()
df5_A <- list()
df5_B <- list()
df5_A2 <- list()
df5_B2 <- list()

# remove burnin
for (i in 1:length(data)) {
  start[[i]] <- round(0.1001*length(data[[i]]$extinction.1))
  end[[i]] <- length(data[[i]]$extinction.1)
}

# get speciation and extinction rates
for (i in 1:length(data)) {
  HiSSE_types[[i]] <- rep(c("0A", "1A", "0B", "1B"), each = length(data[[i]]$extinction.1[start[[i]]:end[[i]]]))
  dat_div[[i]] <- data.frame(dens = c(data[[i]]$speciation.1[start[[i]]:end[[i]]]-data[[i]]$extinction.1[start[[i]]:end[[i]]],
                                      data[[i]]$speciation.2[start[[i]]:end[[i]]]-data[[i]]$extinction.2[start[[i]]:end[[i]]],
                                      data[[i]]$speciation.3[start[[i]]:end[[i]]]-data[[i]]$extinction.3[start[[i]]:end[[i]]],
                                      data[[i]]$speciation.4[start[[i]]:end[[i]]]- data[[i]]$extinction.4[start[[i]]:end[[i]]]),
                             Type = HiSSE_types[[i]])
}

# get diversification rates
for (i in 1:length(data)) {
  # split div based on 0A, 1A, 0B, 1B
  dat_div_split[[i]] <- split(dat_div[[i]], f=dat_div[[i]]$Type)
  A1_0[[i]] <- dat_div_split[[i]][[3]]$dens - dat_div_split[[i]][[1]]$dens
  B1_0[[i]] <- dat_div_split[[i]][[4]]$dens - dat_div_split[[i]][[2]]$dens
}


# keep every 5th step of the MCMC only
for (i in 1:length(data)) {
  A1_0df[[i]] <- data.frame(A1_0[[i]])
  B1_0df[[i]] <- data.frame(B1_0[[i]])
  df5_A[[i]] <- A1_0df[[i]][seq(5, nrow(A1_0df[[i]]), 5), ]
  df5_B[[i]] <- B1_0df[[i]][seq(5, nrow(B1_0df[[i]]), 5), ]
} ## df5_A and df5_B are the diversification for A and for B for every 5th step

# separate extinction-estimated (full) and extinction-constrained (constrained)
full <- c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)
constrained <- c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47)

df5_A_full <- df5_A[full]
df5_B_full <- df5_B[full]

df5_A_c <- df5_A[constrained]
df5_B_c <- df5_B[constrained]

df5_together_full <- list()
df5_together_constrained <- list()

df5_A2_full <- list()
df5_B2_full <- list()
df5_A2_constrained <- list()
df5_B2_constrained <- list()

for (i in 1:length(df5_A_full)) {
  df5_A2_full[[i]] <- data.frame(df5_A_full[[i]])
  df5_B2_full[[i]] <- data.frame(df5_B_full[[i]])
  df5_A2_constrained[[i]] <- data.frame(df5_A_c[[i]])
  df5_B2_constrained[[i]] <- data.frame(df5_B_c[[i]])
  df5_together_full[[i]] <- data.frame(df5_A2_full[[i]], df5_B2_full[[i]])
  df5_together_constrained[[i]] <- data.frame(df5_A2_constrained[[i]], df5_B2_constrained[[i]])
}

# load simmaps to get time spent in each state
files2 <- list.files(pattern="*.character.tree")
data2 <- lapply(files2, read.simmap, format="phylip")

# make empty lists
times_0A <- list()
times_0B <- list()
times_1A <- list()
times_1B <- list()
times_A <- list()
times_B <- list()
times_total <- list()
prop_0A <- list()
prop_1A <- list()
prop_0B <- list()
prop_1B <- list()
propA <- list()
propB <- list()
propA2 <- list()
propB2 <- list()
df5_together_full2 <- list()
df5_together_full3 <- list()
df5_together_constrained2 <- list()
df5_together_constrained3 <- list()

# loop over simmaps to get branches per states and proportions in each
for (i in 1:length(data2)) {
  
  maps <- unlist(data2[[i]]$maps)
  grp <- split(maps, attributes(maps))
  time_0A <- sum(grp$`0`)
  time_1A <- sum(grp$`1`)
  time_0B <- sum(grp$`2`)
  time_1B <- sum(grp$`3`)
  time_A <- sum(time_0A, time_0B)
  time_B <- sum(time_1A, time_1B)
  time_total <- sum(time_A, time_B)
  
  times_0A[[i]] <- time_0A
  times_0B[[i]] <- time_0B
  times_1A[[i]] <- time_1A
  times_1B[[i]] <- time_1B
  times_A[[i]] <- time_A
  times_B[[i]] <- time_B
  times_total[[i]] <- time_total
  
  prop_0A[[i]] <- times_0A[[i]]/times_total[[i]]
  prop_1A[[i]] <- times_1A[[i]]/times_total[[i]]
  prop_0B[[i]] <- times_0B[[i]]/times_total[[i]]
  prop_1B[[i]] <- times_1B[[i]]/times_total[[i]]
  
  propA[[i]] <- prop_0A[[i]] + prop_1A[[i]]
  propB[[i]] <- prop_0B[[i]] + prop_1B[[i]]
}

# convert to data frame
for (i in 1:length(propA)) {
  propA2[[i]] <- data.frame(propA[[i]])
  propB2[[i]] <- data.frame(propB[[i]])
}

# separate extinction-estimated (full) and extinction-constrained (constrained)
full <- c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)
constrained <- c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47)

propA2_full <- propA2[full]
propB2_full <- propB2[full]
propA2_c <- propA2[constrained]
propB2_c <- propB2[constrained]

# join diversification estimates above with simmap outputs
for (i in 1:length(df5_together_full)) {
  df5_together_full2[[i]] <- data.frame(c(df5_together_full[[i]],propA2_full[[i]]))
  df5_together_full3[[i]] <- data.frame(c(df5_together_full2[[i]],propB2_full[[i]]))
  df5_together_constrained2[[i]] <- data.frame(c(df5_together_constrained[[i]],propA2_c[[i]]))
  df5_together_constrained3[[i]] <- data.frame(c(df5_together_constrained2[[i]],propB2_c[[i]]))
  }

for (i in 1:length(df5_together_full)) {
df5_together_full3[[i]]$value <- (df5_together_full3[[i]]$df5_A_full..i..*df5_together_full3[[i]]$propA..i..) + (df5_together_full3[[i]]$df5_B_full..i..*df5_together_full3[[i]]$propB..i..)
df5_together_constrained3[[i]]$value <- (df5_together_constrained3[[i]]$df5_A_c..i..*df5_together_constrained3[[i]]$propA..i..) + (df5_together_constrained3[[i]]$df5_B_c..i..*df5_together_constrained3[[i]]$propB..i..)
}

# get means, medians, and quantiles
means1 <- list()
means2 <- list()
quantiles1 <- list()
quantiles2<- list()
medians1 <- list()
medians2 <- list()

for (i in 1:length(df5_together_full)) {
  means1[[i]] <- mean(df5_together_full3[[i]]$value)
  means2[[i]] <- mean(df5_together_constrained3[[i]]$value)
  quantiles1[[i]] <- quantile(df5_together_full3[[i]]$value)
  quantiles2[[i]] <- quantile(df5_together_constrained3[[i]]$value)
  medians1[[i]] <- median(df5_together_full3[[i]]$value)
  medians2[[i]] <- median(df5_together_constrained3[[i]]$value)
}

q21<-list()
q41<-list()
q22<-list()
q42<-list()
for (i in 1:length(df5_together_full)) {
 q21[[i]] <-  unlist(quantiles1[[i]][[2]])
 q41[[i]] <-  unlist(quantiles1[[i]][[4]])
 q22[[i]] <-  unlist(quantiles2[[i]][[2]])
 q42[[i]] <-  unlist(quantiles2[[i]][[4]])
}

# save weighted diversification rates
full_data <- data.frame(unlist(means1), unlist(medians1), unlist(q21), unlist(q41))
constrained_data <- data.frame(unlist(means2), unlist(medians2),unlist(q22), unlist(q42))
write.csv(full, "full_med.csv")
write.csv(constrained, "constrained_med.csv")

