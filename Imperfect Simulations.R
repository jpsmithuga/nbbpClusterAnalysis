# rm(list=ls())
##'This function simulates a branching process and generates a single transmission chain 
##' n+1 generation given by Z_{n+1} = X^n_1 + X^n_2 + ... + X^n_{Zn} 
##' where X^n_i are IID branching offspring dist 
##' @param gens number of generations to go out to, excluding initial, default to 100 to ensure sufficiently large
##' @param init.size initial size, default to 1 
##' @param offspring distribution to be specified
##' @param ... for passing variable number of inputs depending on choice of offspring
##--------------------------------------------------------------------------------------------------------------
bp <- function(gens = 100, init.size = 1, offspring, ...){  
  Z <- list()
  Z[[1]] <- init.size
  i <- 1 
  while(sum(Z[[i]]) > 0 && i <= gens) { 
    Z[[i+1]] <- offspring(sum(Z[[i]]), ...) 
    i <- i+1 
  } 
  return(Z)
} 

##' _______________________________________________________________________________________________
##' Simulating imperfect observation 
##'      @param true_R       True underlying R value for NB branching process
##'      @param true_k       True underlying k value for NB branching process
##'      @param num_chains   Number of simulated transmission chains in a surveillance system
##'      @param p1           Probability of ascertaining cases by passive surveillance
##'      @param p2           Probability of ascertaining cases by active surveillance
##'      @param prob_cens    Probability that a chain will be censored
##'      @param perc_overlap Proportion of clusters that overlap (i.e. multiple index cases)
##'      - - - - - - - - - - - - -       
##'      @return Output is a list of length 2, each list contains a data frame of cluster sizes, index 
##'              cases, and censored status, for: 
##'                 [[1]] Perfect Surveillance
##'                 [[2]] Imperfect Surveillance
##' _______________________________________________________________________________________________
##' 
true_R <- 0.9
true_k <- 0.25
prob_cens <- perc_overlap <- 0.1

imperfect <- function(true_R, true_k, num_chains = 2000, p1, p2, prob_cens, perc_overlap){
  z <- replicate(num_chains,bp(offspring = rnbinom, mu = true_R, size = true_k)) 
  z.pass <- z; z.act <- z                               
  ## - - - - - - - - - - - - - - - - - - - - -
  ## Imperfect case ascertainment
  ## - - - - - - - - - - - - - - - - - - - - -
  #Passive surveillance
  for (i in 1:length(z.pass)){ 
    for (j in 1:length(z.pass[[i]])){
      for (k in 1:length(z.pass[[i]][[j]])){
        for (l in 1:length(z.pass[[i]][[j]][[k]])){
          if (runif(1) < (1 - p1)){                     
            z.pass[[i]][[j]][[k]] <- NA 
          }}}}}
  #Active case finding (only chains with at least one case observed by passive surveillance)
  ##' Note: Index cases with no secondary cases present in the branching process as a 
  ##' list of 2 (always [1] and [0]). The second value is what we are concerned about having an NA value via passive
  ##' surveillance. In the scenario where the list is [1] and [NA], this would artificially make the chain eligible for
  ##' reevaluation since the cluster is not completely missing. So before reevaluation with p2 we need to 
  ##' assign the first position in all of the lists as NA. That way, if the index case is missing, 
  ##' both values will be missing.
  a <- z.pass
  for (i in 1:length(a)){
    a[[i]][[1]] <- NA
  }
  b <- cbind(a, sapply(a, function(x) all(is.na(unlist(x))))) # Determine if at least one case in the chain has been seen
  for (i in 1:length(z.pass)){
    if (b[i,2] == TRUE){  
      z.act[[i]] <- a[[i]] # Skip if cluster no cases are observed
    } else {              
      for (j in 1:length(z.pass[[i]])){
        for (k in 1:length(z.pass[[i]][[j]])){
          for (l in 1:length(z.pass[[i]][[j]][[k]])){
            if (is.na(z.pass[[i]][[j]][[k]])){      
              if(runif(1) <= (p2)){       #Active probability of being seen by case detection         
                z.act[[i]][[j]][[k]] <- z[[i]][[j]][[k]] #Reassign original value
              } else {z.act[[i]][[j]][[k]] <- z.pass[[i]][[j]][[k]]}
            }}}}}}
  # "Break" chains based on the position of missing cases in the chain
  l <- z.act  #dummy/temp data to not change z.act
  for (i in 1:length(l)){ 
    l[[i]][[1]] <- NULL #remove first position of the nested list so that it eases summing lengths (can't sum based on integer values in imperfect observations)
  }
  #"Break apart" the chains 
  t1 <- lapply(lapply(seq_along(l), function(nm) {split(l[[nm]], cumsum(sapply(l[[nm]], function(x) all(is.na(x)))))}), function(lstA) lapply(lstA,function(x) Filter(function(y) !all(is.na(y)), x)))
  t2 <- rapply(unlist(t1,recursive=FALSE),function(x) x[!is.na(x)], how="replace") #Remove NA values. 
  z.broken <- Filter(length,t2) #remove all with length 0 (missing/unobserved)

  ## - - - - - - - - - - - - - - - - - - - - -
  ## Censoring
  ## - - - - - - - - - - - - - - - - - - - - -
  
  z.cen <- z.broken                          # Initialize the censored list
  for (i in 1:length(z.broken)){             # Iterate through the list
    if (length(z.broken[[i]]) > 1) {           # List must have at least length of two (cant be censored if the index case isn't seen, then it is unobserved as above)
      if(runif(1) <= prob_cens){               # Stochastic process to determine if the nested list will be censored
        if(length(z.broken[[i]]) == 2){
          n <- 2} else {                     # n has to be 2 for chains with only two generations
            n <- sample(2:length(z.broken[[i]]), 1)}   # Randomly determine what list position in the nested list will be the censor threshold
        z.cen[[i]][n:length(z.broken[[i]])] <- NA   # Fill all positions from n to the end of the nested list with NA
      }}}
  out_list <- lapply(z.cen, function(x) {    # Remove all nested list elements that contain NA 
    inds <- sapply(x, function(x) any(is.na(x)))
    if(any(inds)) x[seq_len(which.max(inds) - 1)] else x})
  cens <- numeric(length(out_list))
  true <- numeric(length(out_list))
  for (k in 1:length(out_list)){
    cens[k] <- sum(lengths(out_list[[k]])) # Get cluster size of censored clusters
    true[k] <- sum(lengths(z.broken[[k]])) # Get cluster size of uncensored (but imperfect obs) clusters
  }
  Y_cens <- data.frame(y.cens = cens, censor = ifelse(cens!=true,1,0)) #Create a censoring index (1=censored, 0=uncensored)

  ## - - - - - - - - - - - - - - - - - - - - -
  ## Overlapping clusters
  ## - - - - - - - - - - - - - - - - - - - - -
  #' Determine sampling space - i.e. how many clusters get merged with each iteration
  #'     choose j from Poisson with lambda=2
  lamda <- 2
  a <- rpois(10000000,lamda) # Drawn from a Poisson with lamda=2 for main text
  n <- round(nrow(Y_cens)*perc_overlap)         # Determine the number to be merged
  sample_clust <- sample(a[a > 1], size = n, replace = TRUE) # Randomly choose j (discard 0 and 1)
  names(sample_clust) <- paste0("S", 1:n)
  m <- nrow(Y_cens) - sum(sample_clust) # Determine the number that will not be merged (m)
  non_merge_clust <- rep(1, m)        # Create a vector with replicated 1 based on m
  names(non_merge_clust) <- paste0("N", 1:m)
  combine_clust <- c(sample_clust, non_merge_clust) # Combine sample_clust and non_merge_clust, and then randomly sort the vector
  combine_clust2 <- sample(combine_clust, size = length(combine_clust))
  expand_list <- list(lengths = combine_clust2, values = names(combine_clust2))# Expand the vector
  expand_clust <- inverse.rle(expand_list)
  dat <- data.frame(Y_cens, group = factor(expand_clust, levels = unique(expand_clust)))# Create a data frame with y and expand_clust
  dat$index <- 1 # add the index case number for summing
  dat2 <- aggregate(cbind(dat$y.cens, dat$index, dat$censor), by = list(group = dat$group), FUN = sum)  # Convert dat2 to a matrix, sum the index cases and censoring index, remove the group column
  dat2$group <- NULL   
  y.merged <- as.matrix(dat2); colnames(y.merged) <- c("clust_size","index_cases","censor_status") 
  y.final <- data.frame(y.merged) #just to be safe
  y.final$censor_status <- ifelse(y.final$censor_status >= 1,1,0) # if more than 1 censored clusters merged
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  y.true <- unlist(lapply(z,function(x) sum(unlist(x)))) # Sum true cluster sizes
  Y.true <- data.frame(y.true = y.true,  index = rep(1, times = length(y.true)), censor = rep(0, times = length(y.true)))  # Format data for perfect surveillance (no censoring or subclustering)
  names(Y.true) <- c("clust_size","index_cases","censor_status")
  return(list(Y.true, y.final))
  #return(list(z, z.pass, z.act, z.broken, out_list, cens, true, Y_cens, y.final)) #for validation
}



###############################################################################################
##' Example surveillance system originating with 5000 chains 
##' with transmission assumed negative binomially distributed
##' with R=0.9 and k=0.25
##' Assuming: 
##'     50% of cases are observed (p1=0.50)
##'     25% of unobserved cases are found by active case finding (p2=0.25)
##'     10% of clusters cant be teased apart (overlap; pover=0.10)
##'     10% of clusters are censored (pcens = 0.10)
###############################################################################################
set.seed(05062020)
num_chains <- 5000
R <- 0.9
k <- 0.25
p1 <- 0.90
p2 <- 0.25
pcens <- 0.1
pover <- 0.1

#' The `imperfect()` function returns a list containing both perfect and imperfect surveillance for verification
Y_full <- imperfect(true_R = R, true_k = k, num_chains = num_chains, 
                    p1 = p1, p2 = p2, prob_cens = pcens, perc_overlap = pover)

Y_perfect <- Y_full[[1]]
Y_imperfect <- Y_full[[2]]

#' Note that true clusters of size Y=1 cannot be censored (they are either observed or not), so 
#' imperfect surveillance where Y=1 and cens=1 indicates censoring took place at the second generation




