library(MCMCpack)
library(mvtnorm)
library(dplyr)
mixgauss <- read.table("./mixgauss.dat", header=FALSE)
n <- nrow(mixgauss) #mixgauss is n x 3 shape
p <- 3
ttt <- 0

#' Function get class center
#' 
#' @param X data matrix (n,3)
#' @param Z vector class label
#' @return mu class center
class_center <- function(X, Z){
  Xnew <- as.data.frame(cbind(X, Z))
  
  # Group by class label and calculate the mean for each group
  colnames(Xnew) <- c("x", "y", "z", "class_label")
  
  class_centers <- Xnew %>%
    group_by(class_label) %>%
    summarise(center_x = mean(x),
              center_y = mean(y),
              center_z = mean(z)) %>%
    ungroup() %>%
    arrange(class_label)
  # Return the result as a matrix
  return(as.matrix(class_centers[, c("center_x", "center_y", "center_z")]))
}

#' Gibbs sampler for p, \mu, \Sigma, Z
#' @param X data matrix (n,3)
#' @param Z vector class label
#' @param p vector p1,p2,...pK-1
#' @param Sigma vector of matrix 
#' @param mu matrix, row i is mu_i
#' @return sample of step t
Gibbs <- function(X, Z, p, Sigma, mu, K){
  ni <- 1:K
  ni <- sapply(ni, function(x) sum(Z==x))
  
  
  #sample p
  p.t <- c(rdirichlet(1, 1+ni))
  #sample mu
  mu.t <- matrix(0, nrow=K, ncol=3)
  
  #get class center
  class_centers <- class_center(X, Z)
  for (i in 1:K){
    mu.t[i,] <- c(mvtnorm::rmvnorm(1,
                                   class_centers[i,],
                                   Sigma[[i]]/ni[i]))
  }
  #sample Sigma
  Sigma.t <- vector("list", K)
  for (i in 1:K){
    Xi <- t(X[Z==i, ])
    Si <- tcrossprod(Xi - mu[i,])
    tryCatch({
      # Code that might produce an error or warning
      solve(Si)
    }, error = function(err) {
      print(Si)
      print(Z)
      ttt<<-Si
      print(ttt)
      print(i)
      cat("An error occurred:", conditionMessage(err), "\n")
      # You can also access more information about the error using conditionName(err), conditionCall(err), etc.
    }, warning = function(warn) {
      # Code to handle a warning
      cat("A warning occurred:", conditionMessage(warn), "\n")
    })
    
    
    
    
    
    
    Sigma.t[[i]] <- MCMCpack::riwish(ni[i]+1, Si)
  }
  
  #unnormalized prob
  pZi <- 1:K
  pZi <- sapply(pZi,
                function(i) (dmvnorm(X, mean=mu.t[i, ], sigma=Sigma.t[[i]])*p.t[i])
  )
  #print(pZi)
  Z.t <- 1:nrow(X)
  Z.t <- sapply(Z.t, function(i)
    sample(1:K, size=1, replace=TRUE,
           prob=pZi[i,]))
  
  return(list(Z=Z.t, mu=mu.t, Sigma=Sigma.t, p=p.t))
}

#sampling

sampling_4_fix_K <- function(K, X){
  #we need initial value for  Z_1,\cdots,Z_n, p_1,\cdots,p_{K},
  #\mu_1,\Sigma_1,,\cdots,\mu_K,\Sigma_K,
  #K
  p.0 <- rep(1/K, K) # vector p1,p2,...pK
  Z.0 <- sample(1:K, n, replace=TRUE) # class label Z1...Zn
  mu.0 <- matrix(0, ncol=p, nrow=K) #class center
  Sigma.0 <- vector("list", K)
  for (i in 1:K){
    Sigma.0[[i]] <- diag(p)
  }
  
  samp.size <- 10000
  samples.out <- vector("list", samp.size+1)
  samples.out[[1]] <- list(Z=Z.0, mu=mu.0,
                           Sigma=Sigma.0, p=p.0)
  for (i in 2:(samp.size+1)){
    samples.out[[i]] <- do.call(Gibbs,
                                args=c(list(X=X, K=K), samples.out[[i-1]]))
  }
  
  mu.sample <- array(unlist(lapply(samples.out, function(p) p$mu)),
                     dim = c(K, 3, samp.size+1))
  Z.sample <- do.call(rbind, lapply(samples.out, function(p) p$Z))
  print(paste('size:', dim(Z.sample)))
  Sigma.sample <- lapply(samples.out, function(p) p$Sigma)
  
  #point estimate
  #use last 5000 sample
  mu.hat <- rowMeans(mu.sample[ , , 5000:samp.size+1], dims=2)
  
  #zhat we use round
  Z.hat <- round(colMeans(Z.sample[5000:samp.size+1, ]))
  
  Sigma.hat <- vector("list", K)
  for (i in 1:K){
    Sigmai <- lapply(Sigma.sample, function(p) p[[i]])
    Sigmai.hat <- array(unlist(Sigmai),
                        dim=c(3,3, samp.size+1))
    Sigma.hat[[i]] <- rowMeans(Sigmai.hat[ , , 5000:samp.size+1], dims=2)
  }
  
  #BIC
  #-2log(y|theta)+klog(n)
  #likelihood
  #&\prod_{i=1}^np(X_i|Z_i, \boldsymbol{\mu}, \mathbf{\Sigma}) \\
  #=&\prod_{t=1}^K \prod_{i:Z_i=t} f(X_i|\mu_t,\Sigma_t)
  likelihood <- 1:K
  likelihood <- -2*sum(log(sapply(likelihood,
                                  function(t) prod(dmvnorm(X[Z.hat==t,],
                                                           mean=mu.hat[t,],
                                                           sigma=Sigma.hat[[t]])))))
  bic <- likelihood + log(nrow(X))*(
    K*3+K*6+K-1 #mu #sigma #p
  )
  
  return (list(bic=bic,
               mu.hat=mu.hat,
               Z.hat=Z.hat,
               Sigma.hat=Sigma.hat))
}

print(sampling_4_fix_K(X=mixgauss, K=2))
