##### Diversification Rate from HiSSE #####

## get speciation and extinction rates after transforming HiSSE models,
## input rates into equations from Caetano et al 2018

## load data
hisse_model <- readRDS(file="1_full_hisse.rds")

# number of tips with trait state 0:
n0 <- dim(hisse_model$data)[1] - sum(as.numeric(hisse_model$data[,2]))
# number of tips with trait state 1:
n1 <- sum(as.numeric(hisse_model$data[,2]))

## transform turnover and extinction_fraction to speciation and extinction
# speciation = turnover / (1 + extinction_fraction)
# extinction = (turnover*extinction_fraction) / (1 + extinction_fraction)

# speciation0A:
sp0A <- hisse_model$solution[1]/(1 + hisse_model$solution[3])
# extinction0A:
e0A <- (hisse_model$solution[1]*hisse_model$solution[3])/(1 + hisse_model$solution[3])
# speciation0B:
sp0B <- hisse_model$solution[13]/(1 + hisse_model$solution[15])
# extinction0B:
e0B <- (hisse_model$solution[13]*hisse_model$solution[15])/(1 + hisse_model$solution[15])

# speciation1A:
sp1A <- hisse_model$solution[2]/(1 + hisse_model$solution[4])
# extinction1A:
e1A <- (hisse_model$solution[2]*hisse_model$solution[4])/(1 + hisse_model$solution[4])
# speciation1B:
sp1B <- hisse_model$solution[14]/(1 + hisse_model$solution[16])
# extinction1B:
e1B <- (hisse_model$solution[14]*hisse_model$solution[16])/(1 + hisse_model$solution[16])

# div 0A:
div0A <- sp0A-e0A
# div 0B:
div0B <- sp0B-e0B
# div 1A:
div1A <- sp1A-e1A
# div 1B:
div1B <- sp1B-e1B
