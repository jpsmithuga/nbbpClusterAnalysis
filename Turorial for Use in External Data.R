# rm(list=ls())
##' This code will explain how to use these methods with external surveillance cluster data
##' 
################################################################################################
## Preparing the data
################################################################################################

##' First, the data need to be prepared with three columns:
##'     1: cluster size (with each row a unique cluster)
##'     2: number of subclusters/index cases
##'     3: censorship status
##'
##' Typically, the second and third columns are 1 and 0, respectively

##' Each row in the dataset needs to be an individual cluster, so the
##' length of the dataset is the number of clusters in the surveillance
##' system. However, often data come in the form of a table, i.e.:
##' ----------------------------
##' Cluster size | Frequency
##' ----------------------------
##'            1 | 4003
##'            2 | 362
##'            3 | 174
##'            4 | 92
##'            5 | 75
##'            6 | 41
##'            7 | 44
##'            8 | 19
##'            .....
##'            52 | 1
##'            68 | 1
##'            77 | 1
##' 
##' The below code will help you convert table data of cluster sizes to the
##' format appropriate for this code 

##' First, we will create example table data for demonstration purposes using the branching 
##' process function to generate a surveillance system

bp <- function(gens = 100, init.size = 1, offspring, ...){  
  Z <- list() #initiate the list
  Z[[1]] <- init.size #set the first position of the list as the number of index cases
  i <- 1 
  while(sum(Z[[i]]) > 0 && i <= gens) { 
    Z[[i+1]] <- offspring(sum(Z[[i]]), ...) 
    i <- i+1 
  } 
  return(Z)
} 
##' We can now generate our example table data. Here we will simulate a 
##' surveillance system originating with 1000 chains, with an offspring
##' distributed with mean R = 0.50 and dispersion k = 0.15 

num_chains <- 1000
R <- 0.50
k <- 0.15

set.seed(05062020)
table_data <- data.frame(table(unlist(lapply(replicate(num_chains, bp(offspring = rnbinom, mu = R, size = k)),function(x) sum(unlist(x))))))
names(table_data) <- c("cluster_size","frequency") # rename columns
table_data[ , "cluster_size"] <- as.numeric(as.character(table_data[ , "cluster_size"])) # convert to numeric

##' In these data, the first column is the cluster size, the second column is the frequency. 
##' Ensure your data are in this format.
##' To convert it to the proper long format (each individual cluster is a row) for use in these methods, 
##' run the following step on your table data. This assumes no subclustering or censoring.

Y_data <- data.frame(clust_sizes = table_data[rep(seq_len(nrow(table_data)), table_data[, "frequency"]), ][ , 1],
                     index_cases = rep(1, times = sum(table_data[, "frequency"])), 
                     cens_status = rep(0, times = sum(table_data[, "frequency"])))

##' If you want to assume some definition for censoring, i.e., clusters larger than some threshold are censored, 
##' you can easily do this:
##' 
# cens_threshold <- 10
# Y_data[Y_data[ , "clust_sizes"] >= cens_threshold, "cens_status"] <- 1 # all clusters Y >= 10

################################################################################################
## Using the likelihood and parameter estimation functions
################################################################################################
##' Now that the data are prepared and in the proper format, parameter estimation is quite easy
##' Pull in the functions from the file, "Likelihood and Parameter Estimation Functions.R"
##' Note: the functions are also pasted below for convenience, in case the additional file is unavailable

# source("~Insert_Filepath/Likelihood and Parameter Estimation Functions.R")

##' If not yet defined, define search grid, which will be an R x k matrix with each cell containing 
##' the likelihood of each R/k combination

resolution <- 0.01 # set resolution (increased resolution increases computer time needed)

# R range
R.min <- 0.01
R.max <- 1.00
Rrange <- seq(R.min, R.max, by = resolution)

# k range
k.min <- 0.01
k.max <- 1.00
krange <- seq(k.min, k.max, by = resolution)

# Calculate grid of likelihoods 
Y_surflikes <- surflike(data = Y_data, Rrange = Rrange, krange = krange)

# find maximum in the x,y grid
Y_maxlikes <- Y_surflikes == max(Y_surflikes)

# estimate parameters
Y_MLE <- calc_profile(ls = Y_surflikes, ls_max = Y_maxlikes, Rrange = Rrange, krange = krange, conf.interval = 95)
Y_MLE

################################################################################################
## Creating a figure
################################################################################################

# Set values to calculate 90 and 95% confidence regions
CI95 <- qchisq(0.95, df = 1) / 2 
CI90 <- qchisq(0.90, df = 1) / 2 

# Calculate contour lines for 90 and 95% confidence regions
ctlns_95 <- contourLines(x = Rrange, y = log10(krange), as.matrix(Y_surflikes - max(Y_surflikes)), levels = c(-CI95, CI95+0.001))
ctlns_90 <- contourLines(x = Rrange, y = log10(krange), as.matrix(Y_surflikes - max(Y_surflikes)), levels = c(-CI90, CI90+0.001))

# Set some general figure parameters, i.e. axes, colors, etc
xtick <- seq(0, 1.4, 0.1)
ytick <- c(0.01, 0.05, 0.1, 0.15, 0.3, 0.5, 0.75, 1.0)

# Make figure
# pdf(paste0("~filepath/MLE_90_and_95_Estimates_",Sys.Date(),".pdf"), width=10, height=8) # if you want to save PDF, uncomment this and "dev.off()" below
filled.contour(x = Rrange, y = log10(krange), as.matrix(Y_surflikes - max(Y_surflikes)), col='white', levels = seq(-3,0,0.5), 
               xlab = 'Reproduction number, R', ylab = "Dispersion parameter, k (log scale)", 
               xlim=c(0.01, 1),
               ylim=log10(c(0.01,max(krange))),
               plot.axes = {
                 axis(1, at = xtick, label = xtick)
                 axis(2, at = log10(ytick), label = ytick)
                 
                 points(Y_MLE["R","point_est"], log10(Y_MLE["k","point_est"]), pch = 19, cex = 0.5, col="black")
                 lines(ctlns_95[[1]][[2]], ctlns_95[[1]][[3]], lty = 1, col = "black") 
                 lines(ctlns_90[[1]][[2]], ctlns_90[[1]][[3]], lty = 3, col = "black")
               }
)
legend('topleft', c("95% CI", "90% CI"), lty = c(1,3), lwd = 2, col = "black", bty='n')
# dev.off()




##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' Functions to estimate R and k from transmission cluster data (pasted for convenience)
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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


