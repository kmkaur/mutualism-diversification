##### Run Sister Contrasts #####

##### Load packages #####
require(ape)
require(phytools)
require(phangorn)
require(dplyr)

##### Load functions for sister contrasts #####
## FUNCTION 1: drop sister tips that share the same trait value
same_nodes <- function(phylogeny, trait, label_tip){
    # get tip numbers
    ntips <- Ntip(phylogeny)
    
    # matrix to fill in with parent nodes (node paths)
    parent_nodes <- matrix(NA_integer_, nrow = ntips, ncol = ntips-1,
    dimnames = list(phylogeny$tip.label,
    seq(ntips+1, (2*ntips)-1, by = 1)))
    # retrieve paths
    full_paths <- path_vec(1:ntips, phylogeny$edge, ntips+1)
    
    # fill in matrix
    for(i in 1:ntips){
        parent_nodes[i, full_paths[[i]]-ntips] =
        trait[which(phylogeny$tip.label == label_tip[i])]
    }
    
    # nodes that will be dropped
    drop_nodes <- apply(parent_nodes, 2, function(x) length(unique(na.exclude(x))) )
    
    return(as.numeric(names(which(drop_nodes == 1))))
}

## FUNCTION 2: retrieve the node numbers for each tip from tip to root (path)
path_retrieve <- function(tip_number, edge, root_number){
    
    #empty list
    final_path <- c()
    
    #get path
    while(tip_number != root_number){
        final_path <- c(which(edge[ , 2] == tip_number), final_path)
        tip_number <- edge[final_path[1], 1]
    }
    
    final_path <- edge[final_path , 1]
    
    return(rev(final_path))
}

# tips as vector
path_vec <- Vectorize(FUN = path_retrieve, vectorize.args = "tip_number")

## NOTE: The following functions are from Käfer and Mousset, 2014 Dryad Repository, with some edits.
# https://doi.org/10.5061/dryad.jd8vg
# Input for the Käfer and Mousset, 2014 paper requires a data table with two columns,
# the first column is 'd' or the number in the derived state (1 in our case) and
# the second column is "m" or the total number in the clade (0's + 1s).
# Each row represents a sister pair.

## Conditional probability of having S>=s given M=m and K=k
PSmk <- Vectorize(function(s, m, k) {
    return( exp( lchoose(m-s, k-1) - lchoose(m-3, k-1) ) )
})

## Conditional Expectation of L (length of the root edge with p desendants) given M=m
ELmp <- Vectorize(function(m,p) {
    a <- (m-p-1)/(m-2)
    b <- (p-1)/2/(m-2)
    if (a>0) {
        s <- 3:(m-p+1) #not sure what s is, 3:total-descendants+1 ?
        psmk <- PSmk(s, m, p)
        a <- a*(0.5 + sum(psmk/s))
    }
    return(a+b)
})

## Conditional probability of D (size of the derived clade) given M=m and K in {d,m-d}
PDmp <- Vectorize( function(d,m) {
    if (d != m-d) {
        a <- ELmp(m,d)
        b <- ELmp(m,m-d)
        return( a/(a+b) )
    } else {
        return(1)
    }
})


## Käfer and Mousset Test for Sister Clade Comparison
## iter  : number of resamplings
scc.test <- function(dataset, iter=10000, alternative="two.sided") {
    
    # Select pairs of clades with m > 2
    # This removes all 1:1 contrasts
    # But can be omitted and still works
    dataset <- dataset[dataset[["m"]]>2,]
    
    # resampling function to generate n simulated d/m values
    resample <- function(y, n=1) {
        f <- function() {
            p <- runif(dim(y)[1]) < y[,2]
            return(mean(c(y[p,1],1-y[!p,1])))
        }
        return(replicate(n=n, expr=f()))
    }
    
    # Generate the matrix for the resampling function
    dpmatrix <- matrix(data=c(dataset[["d"]]/dataset[["m"]],
    PDmp(d=dataset[["d"]], m=dataset[["m"]])),
    ncol=2, byrow=FALSE)
    
    # Observed and simulated mean d/m values
    dmobs <- mean(dpmatrix[,1])
    dsim <- NULL
    nqueue <- NULL
    
    # Iterate the resampling method
    for (it in c(iter[1],diff(iter))) {
        # Sample simulated d/m
        dsim <- c(dsim, resample(dpmatrix, n=it))
        ntotal <- length(dsim)
        # Compute the length of the queue of the distribution
        if (alternative == "two.sided") {
            nupqueue <- sum(dsim >= dmobs)
            nlowqueue <- sum(dsim <= dmobs)
            nqueue <- min(nupqueue, nlowqueue)
        } else if (alternative == "greater") {
            nqueue <- sum(dsim >= dmobs)
        } else if (alternative == "less") {
            nqueue <- sum(dsim <= dmobs)
        }
    }
    # compute the p-value of the test
    if (alternative == "two.sided") {
        p.value <- 2*nqueue/length(dsim)
    } else {
        p.value <- nqueue/length(dsim)
    }
    # return the result of the test
    scc <- list(statistic=dmobs, p.value=p.value, alternative=alternative,  iterations=iter )
    return(scc)
}

##### Load Data #####
tree <- read.tree(file="tree_file")
data <- read.csv(file="traits_file")

#match data and order
to_remove <- setdiff(data$species,tree$tip.label)
tree <- drop.tip(tree, to_remove)
data <- data[!data$species %in% to_remove, ]
data <- data[match(tree$tip.label, data$species),]

##### Contrasts#####
# run conrast function
dtest <- same_nodes(tree, data$data, data$species)
ltest <- lapply(dtest, function(x) extract.clade(tree, x)$tip.label)
rtest <- unique(unlist(sapply(ltest, '[', -1, simplify = TRUE)))
ptest <- drop.tip(phy = tree, tip = rtest)

# save files
saveRDS(dtest, file="1_dtest.rds")
saveRDS(ltest, file="1_ltest.rds")
saveRDS(rtest, file="1_rtest.rds")
saveRDS(ptest, file="1_ptest.rds")

##### Run sister test #####

# load data again with updated format
named_with_trait <- as.data.frame(cbind(paste(data$data,"-",data$species), paste(data$species)))

# retrieve the sister pairs
descendants_list <-lapply(1:ptest$Nnode+Ntip(ptest),function(n,t) phangorn::Descendants(t,n)[[1]],t=ptest)
nodes_list <- c(1:ptest_ofs_1090$Nnode+Ntip(ptest))[which(sapply(descendants_list,length)==2)]
sisters_list <- t(sapply(nodes_list,function(n,t) ptest$tip.label[phangorn::Descendants(t,n)[[1]]],t=ptest))
sisters_data_frames <- as.data.frame(sisters_list)
colnames(sisters_data_frames) <- c("tip1", "tip2")
numbers_contrasts <- length(sisters_data_frames$tip1) #need this output summary
colnames(named_with_trait) <- c("tip_new", "tip_old")
test_replace <- sisters_data_frames
test_replace <- as.data.frame(lapply(test_replace , function(x) named_with_trait$tip_new[match(x, named_with_trait$tip_old)]))
colnames(test_replace) <- c("tip0", "tip1")
test_replace$ordering <- ifelse(grepl('0 -', test_replace1090$tip0), 'YES', 'NO')
test_replace[test_replace0$ordering == "NO", c("tip0", "tip1")] <- test_replace[test_replace_sam$ordering == "NO", c("tip1", "tip0")]
test_replace <- test_replace[!grepl("0 - ", test_replace$tip1),]
test_replace <- test_replace[!grepl("1 - ", test_replace$tip0),]
test_replace$ordering <- NULL
test_replace$tip0 <- gsub("0 - ", "", test_replace$tip0)
test_replace$tip1 <- gsub("1 - ", "", test_replace$tip1)

# retrieve the number of nodes represented by the sister pairs
t_remain <- lapply(ltest , `[[`, 1)
t_remain2 <- sapply(ltest, '[', -1, simplify = TRUE)
t_remain3 <- lapply(t_remain2, function(x) length(x))
t_and_r <- data.frame(unlist(t_remain), unlist(t_remain3))
colnames(t_and_r) <- c("tip", "number")
t_and_r<- t_and_r %>% arrange(desc(number))  # desc orders from largest to smallest
tr_unique <- t_and_r  %>% distinct(tip, .keep_all = TRUE)
p_i_d <- tr_unique[tr_unique$tip %in% ptest$tip.label, ]
p_i_d$number <- p_i_d$number + 1
a_t <- setdiff(ptest_ofs$tip.label, p_i_d$tip)
ones_vector <- replicate(length(a_t), 1)
add_this <- data.frame(a_t, ones_vector)
colnames(add_this) <- c("tip", "number")
pid_all <- rbind(p_i_d, add_this)

# make the sister tables with the numbers and run the sister test
colnames(pid_all) <- c("tip0", "number")
sister_0_numbers <- merge(test_replace, pid_all[, c("tip0", "number")], by="tip0")
colnames(pid_all) <- c("tip1", "number")
sister_1_numbers <- merge(test_replace, pid_all[, c("tip1", "number")], by="tip1")
sister_0_numbers$tip0 <- sapply(sister_0_numbers$tip0, function(x) gsub("pEMP_", "", x))
sister_1_numbers$tip0 <- sapply(sister_1_numbers$tip0, function(x) gsub("pEMP_", "", x))
sister_0_numbers <- sister_0_numbers[gtools::mixedorder(sister_0_numbers$tip0) , ]
sister_1_numbers <- sister_1_numbers[gtools::mixedorder(sister_1_numbers$tip0) , ]
sisters_both_numbers <- data.frame(sister_0_numbers$number, sister_1_numbers$number)
sisters_both_numbers$difference <- sisters_both_numbers$sister_0_numbers.number - sisters_both_numbers$sister_1_numbers.number
sisters_both_numbers$total <- sisters_both_numbers$sister_0_numbers.number + sisters_both_numbers$sister_1_numbers.number
tally_differences <- data.frame(t(data.frame(table(sign(sisters_both_numbers$difference)))))
sister_test <- ape::diversity.contrast.test(sisters_both_numbers, method="difference") #need this output summary
wt_sister_test <- wilcox.test(sisters_both_numbers$sister_0_numbers.number, sisters_both_numbers$sister_1_numbers.number, paired=TRUE, exact=FALSE)

# run the resampling on the sister pairs
for_resample <- cbind(sisters_both_numbers$sister_1_numbers.number, sisters_both_numbers$total)
colnames(for_resample) <- c("d", "m") #where d is number '1' or derived, and m is total in both sisters
for_resample <- as.data.frame(for_resample)
scc.test(for_resample, iter = 10000, alternative = "two.sided")

# save files
saveRDS(test_replace_sam, file="sister_contrasts.rds")
saveRDS(sisters_both_numbers, file="sister_numbers.rds")

