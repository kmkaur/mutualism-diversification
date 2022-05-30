##### corHMM & Simmap #####
# load packages
require(corHMM)
require(ape)
require(phytools)
require(hisse) # 1.9.18 is the newest package version
require(diversitree)
require(phangorn)
require(dplyr)

# load hisse
hisse_1 <- readRDS(file="1_full_hisse.rds")

# extract tree and data (in phytools format)
phy <- hisse_1$phy
data1 <- as.data.frame(hisse_1$data)
data <-setNames(data1$data,data1$species)

## run corHMM
# use the transition rate matrix from hisse in corHMM for stochastic mapping
cm1 <- corHMM(phy = phy, data = data1, rate.cat=2, model="ARD", rate.mat = hisse_1$trans.matrix)
saveRDS(cm1, file="1_corhmm.rds")

# save corHMM output
model1 <- cm1$solution
saveRDS(model1, file="1_corhmm_solution.rds")

## make Simmap
# plot corHMM output into simmap (1,000 simulations)
cm2 <- makeSimmap(tree = cm1$phy, data = cm1$data, rate.cat=2, nSim=1000, model = model1)
saveRDS(cm2, file = "1_simmaps.rds")

## extract time spent in each state according to simmaps

# sum frequency of branches per states 0 and 1 and A and B
# empty lists
times_0A <- list()
times_0B <- list()
times_1A <- list()
times_1B <- list()
times_0 <- list()
times_1 <- list()
times_total <- list()

# loop over 1,000 simmaps
for (i in 1:length(cm2)) {
  
  maps <- unlist(cm2[[i]]$maps)
  grp <- split(maps, attributes(maps))
  time_0A <- sum(grp$`1`)
  time_0B <- sum(grp$`2`)
  time_1A <- sum(grp$`3`)
  time_1B <- sum(grp$`4`)
  time_0 <- sum(time_0A, time_0B)
  time_1 <- sum(time_1A, time_1B)
  time_total <- sum(maps)
  
  times_0A[[i]] <- time_0A
  times_0B[[i]] <- time_0B
  times_1A[[i]] <- time_1A
  times_1B[[i]] <- time_1B
  times_0[[i]] <- time_0
  times_1[[i]] <- time_1
  times_total[[i]] <- time_total
  
}

# mean across simmaps
final_0A <- mean(unlist(times_0A))
final_0B <- mean(unlist(times_0B))
final_1A <- mean(unlist(times_1A))
final_1B <- mean(unlist(times_1B))
final_0 <- mean(unlist(times_0))
final_1 <- mean(unlist(times_1))
final_total <- mean(unlist(times_total))

final_0A
final_0B
final_1A
final_1B

