##' _______________________________________________________________________________________________
##' Branching process 
##'    General function to generate data from branching process with variable offspring assumptions 
##' _______________________________________________________________________________________________
##' 
bp <- function(gens=20, init.size=1, offspring, ...){  
  Z <- list() #initiate the list
  Z[[1]] <- init.size #set the first position of the list as the number of index cases
  i <- 1 
  while(sum(Z[[i]]) > 0 && i <= gens) { 
    Z[[i+1]] <- offspring(sum(Z[[i]]), ...) 
    i <- i+1 
  } 
  return(Z)
} 
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' Functions to estimate R and k from transmission cluster data 
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' _______________________________________________________________________________________________
##' Likelihood Function
##' For use in parameter estimation 
##'      @param Y 3-column data frame or matrix containing 
##'                  [,1] Custer Size
##'                  [,2] Index Cases
##'                  [,3] Censored status
##'                  Where each row is a unique cluster 
##'      @param R NB R value to be plugged into the likelihood
##'      @param k NB k value to be plugged into the likelihood
##'      - - - - - - - - - - - - -       
##'      @return Sum of the log-likelihoods 
##' _______________________________________________________________________________________________

likelihood <- function(Y,R,k) {
  p_function <- function(y,n){         
    exp(log(n)-log(y)+lgamma(k*y+y-n)-(lgamma(k*y)+lgamma(y-n+1))+(y-n)*log(R/k)-(k*y+y-n)*log(1+R/k))
  }
  ya <- Y[Y[,3] == 0,] # Uncensored
  yb <- Y[Y[,3] == 1,] # Censored
  
  liks_a <- log(p_function(ya[,1], ya[,2])) 
  liks_b <- numeric(nrow(yb))              
  
  if(nrow(yb) > 0){                          
    for (i in 1:nrow(yb)){
      y <- yb[i, 1]
      n <- yb[i, 2]
      if (y == 1){                           
        liks_b[i] <- 0                     
      } else{
        liks_b[i] <- log(max(10^-300, 1 - sum(p_function(1:(y-1),n)), na.rm=TRUE)) 
      }}}
  sumliks <- sum(liks_a, liks_b)
  return(sumliks)
  #return(cbind(liks_a,liks_b))
}

##' _______________________________________________________________________________________________
##' Surface/profile likelihood function
##' Calculates likelihoods over a range of R and k values
##'      @param data 3-column data frame or matrix containing 
##'                  [,1] Custer Size
##'                  [,2] Index Cases
##'                  [,3] Censored status
##'      @param Rrange Range of R values
##'      @param krange Range of k values
##'      - - - - - - - - - - - - -       
##'      @return Rrange by krange Matrix with likelihoods
##' _______________________________________________________________________________________________

surflike <- function(data, Rrange, krange){
  likesurf <- matrix(NA, nrow = length(Rrange),length(krange))
  for(i in 1:length(Rrange)){
    for(j in 1:length(krange)){
      likesurf[i,j] <- likelihood(data, Rrange[i],krange[j])
    }
  }
  return(likesurf)
}

##' _______________________________________________________________________________________________
##' Parameter Estimation
##' Estimates MLE and confidence interval for R and k
##'      @param ls likelihood surface data
##'      @param ls_max logical likelihood surface data identifying max  (ls_max <- ls==max(ls))
##'      @param conf.interval Desired confidence interval (as decimal, i.e. 0.95)
##'      - - - - - - - - - - - - -       
##'      @return Point, lower, and upper bound estimates for R and k 
##' _______________________________________________________________________________________________

calc_profile <- function(ls, ls_max, Rrange, krange, conf.interval){
  chiV <- qchisq(conf.interval/100, df = 1) / 2
  prfk <- apply(ls,2,function(x){max(x)})
  prfk2 <- krange[prfk - max(prfk) >- chiV]
  prfR <- apply(ls,1,function(x){max(x)})
  prfR2 <- Rrange[prfR - max(prfR) >- chiV]
  
  output <- rbind(cbind(Rrange[sum(seq(1, length(Rrange)) %*% ls_max)], min(prfR2),max(prfR2)),
                  cbind(krange[sum(ls_max %*% seq(1, length(krange)))], min(prfk2),max(prfk2)))
  colnames(output) <- c("point_est","lower_ci","upper_ci")
  rownames(output) <- c("R","k")
  return(output)
}


######################################################################################################
## Example
######################################################################################################
#' Simulate a surveillance system with 10000 chains with underlying R=0.5 and k=0.15, 
#' assuming perfect surveillance
num_chains <- 10000
R <- 0.50
k <- 0.15

######################################################
## Simulate surveillance data
######################################################
set.seed(05062020)
# Individual-level data (Z values) - i.e. full distribution of individual secondary cases 
Z_values <- replicate(num_chains, bp(offspring = rnbinom, mu = R, size = k))

# Cluster data (Y values) - vector of integers representing final chain size
Y_values <- unlist(lapply(Z_values,function(x) sum(unlist(x))))

######################################################
## Maximum Likelihood Estimation
######################################################

#####################
### Individual data
### For individual level data, use classical methods
library(MASS)
Z_MLE <- fitdistr(unlist(lapply(Z_values, function(x) x[-1])),"Negative Binomial")

#####################
### Cluster data
### For cluster level data, employ the above custom functions for cluster-based MLE
  
# prep data in proper format (assuming perfect surveillance so all n=1 and no censoring)
Y_data <- data.frame(clust_size = Y_values,
                     index_cases = 1,
                     cens_status = 0)

# Define search grid
resolution <- 0.01 # set resolution (increased resolution increases computer time needed)
  # R range
R.min <- 0.01
R.max <- 1
Rrange <- seq(R.min, R.max, by = resolution)
  # k range
k.min <- 0.01
k.max <- 1
krange <- seq(k.min, k.max, by = resolution)

  # Calculate grid of likelihoods 
Y_surflikes <- surflike(data = Y_data, Rrange = Rrange, krange = krange)
  # find max in the x,y grid
Y_maxlikes <- Y_surflikes == max(Y_surflikes)
  # estimate parameters
Y_MLE <- calc_profile(ls = Y_surflikes, ls_max = Y_maxlikes, Rrange = Rrange, krange = krange, conf.interval = 95)




