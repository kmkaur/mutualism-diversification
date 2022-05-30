##### HiSSE Models #####

##### Load packages #####
require(ape)
require(phytools)
require(hisse) # 1.9.18 is the newest package version
require(diversitree)
require(phangorn)
require(dplyr)

##### HiSSE #####
# total tips in the tree
total_tree <- length(tree$tip.label)
# total species in group
# based on estimate
total_sp <- # fill in estimate here

#global:
samp1 <- c(total_tree/total_sp)

#### HiSSE ###
##Is mutualism associated with diversification?

# make tree ultrametric
pruned1 <- phangorn::nnls.tree(cophenetic(tree),tree,rooted=TRUE)

# convert data to hisse format
taxa2 <- as.matrix(data)

## HiSSE models

## BiSSE like HiSSE model
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits = 0)
bisse_like_hisse <- hisse(pruned1, taxa2, hidden.states=FALSE, f=samp1, turnover=c(1,2),
eps=c(1,2), trans.rate=trans.rates.bisse)

## null BiSSE model
null_bisse_like_hisse <- hisse(pruned1, taxa2, hidden.states=FALSE, f=samp1, turnover=c(1,1),
eps=c(1,1), trans.rate=trans.rates.bisse)

## CID-2 HiSSE model
trans.rates.hisse.1 <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
null_two_hisse <- hisse(pruned1, taxa2,f=samp1, hidden.states=TRUE, turnover=c(1,1,2,2),
eps=c(1,1,2,2), trans.rate=trans.rates.hisse.1)

## CID-4 HiSSE model
trans.rates.hisse.2 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
null_four_hisse <- hisse(pruned1, taxa2, f=samp1, hidden.states = TRUE, turnover=c(1, 1, 2, 2, 3, 3,4, 4),
eps=c(1, 1, 2, 2, 3, 3, 4, 4), trans.rate = trans.rates.hisse.2)

## full HiSSE model
trans.rates.hisse.3 <- TransMatMakerHiSSE(hidden.traits=1)
full_hisse <- hisse(pruned1, taxa2, f=samp1, hidden.states = TRUE,turnover=c(1,2,3,4),
eps=c(1,2,3,4), trans.rate = trans.rates.hisse.3)

## HiSSE with one hidden state only
hisse_1hidden <- hisse(pruned1, taxa2, f=samp1, hidden.states=TRUE, turnover=c(1,2,1,2),
eps=c(1,2,1,2), trans.rate=trans.rates.hisse.3)

hisse_1hidden2 <- hisse(pruned1, taxa2, f=samp1, hidden.states=TRUE, turnover=c(1,2,0,3),
eps=c(1,2,0,3), trans.rate=trans.rates.hisse.3)

# Save outputs
saveRDS(bisse_like_hisse, file = "1_bisse_like_hisse.rds")
saveRDS(null_bisse_like_hisse , file = "1_null_bisse_like_hisse.rds")
saveRDS(null_two_hisse, file = "1_null_two_hisse.rds")
saveRDS(null_four_hisse, file = "1_null_four_hisse.rds")
saveRDS(full_hisse, file = "1_full_hisse.rds")
saveRDS(hisse_1hidden, file = "1_hisse_1hidden.rds")
saveRDS(hisse_1hidden2, file = "1_hisse_1hidden2.rds")



