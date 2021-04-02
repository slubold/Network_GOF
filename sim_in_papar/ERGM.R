library(ergm)
library(sna)
library(RMTstat)

n = 100 # graph size
numSim = 100 # number of simulations to estimate pHat
numBoot = 500
null = F
theta = c(-1, -0.5)

# Generate the data 
TestNull <- function(theta, n, numBoot, numSim, null = F) {
  
  t_Phat <- rep(0, numSim)
  t_ERGM <- rep(0, numSim)
  
  for(index in 1:numSim){

    g = simulate(network(n, directed = F) ~ edges + kstar(2), coef = theta, directed = F)
    
    if (null == T){
      MLE.fit <- ergm(g ~ edges + kstar(2)) # fit model
    }
    
    if (null == F){
      MLE.fit <- ergm(g ~ edges) # fit model
    }
    
    G = unname(as.matrix.network.adjacency(g))
    
    # Now simulate from the distribution
    
    if(index == 1){
      
      sumA = matrix(rep(0, n^2), ncol = n)
      while(sum(which(sumA[upper.tri(sumA)] == 0)) > 0){
        
        for(indexSim in 1:numBoot){
          g.sim <- simulate(MLE.fit, nsim = 1, burnin = 1000, interval = 1000)
          A = unname(as.matrix.network.adjacency(g.sim))
          sumA = sumA + A
        }
      }
      
      Phat = sumA / numBoot # estimate the probability matrix.
      sumA = matrix(rep(0, n^2), ncol = n)
      
      # We now also simulate the true probability matrix 
      while(sum(which(sumA[upper.tri(sumA)] == 0)) > 0){
        for(indexSim in 1:numBoot){
          g = simulate(network(n, directed = F) ~ edges + kstar(2), coef = theta, directed = F)
          A = unname(as.matrix.network.adjacency(g))
          sumA = sumA + A
        }
      }
      P = sumA / numBoot
      
    }

    # Do bootstrapping here. 
    A = (G - P) / sqrt((n - 1) * P * (1 - P))
    diag(A) = rep(0, n)
    hatA = (G - Phat) / sqrt((n - 1) * Phat * (1 - Phat))
    diag(hatA) = rep(0, n)
    
    ## First method, using Phat.
    #t_Phat[index] = Bootstrap_Phat(hatA, n, numBoot, Phat)
    t_ERGM[index] = Bootstrap_ERGM(hatA, n, numBoot, Phat, MLE.fit)
    difference = n^(2/3) * (eigen(A)$values[1] - eigen(hatA)$values[1])
  }
  return(c(t_Phat, t_ERGM, MLE.fit$coef, difference))
}

values = TestNull(theta, n, numBoot, numSim, null)
write.csv(values, "ERGMS_Null_TRUE.csv")

Bootstrap_Phat = function(hatA, n, numBoot, Phat){
  
  eig <- c(eigen(hatA)$values[1], eigen(hatA)$values[n])
  eig.btsp <- matrix(rep(0, numBoot*2), ncol = numBoot)
  
  
  for (i in 1:numBoot) {
    G.btsp = matrix(rep(0, n^2), ncol = n) 
    G.btsp[upper.tri(G.btsp)] = runif(choose(n, 2)) < Phat[upper.tri(Phat)]
    G.btsp = G.btsp + t(G.btsp)
    
    A.btsp = (G.btsp - Phat) / sqrt( (n - 1) * Phat * (1 - Phat)) 
    diag(A.btsp) = 0
    
    eig.btsp[1, i] = eigen(A.btsp)$values[1]
    eig.btsp[2, i] = eigen(A.btsp)$values[n]  
  }
  
  mu.max = mean(eig.btsp[1,])
  mu.min = mean(eig.btsp[2,])
  sd.max = sd(eig.btsp[1,])
  sd.min = sd(eig.btsp[2,])
  t = -1.2065335745820 + sqrt(1.607781034581) * max( (eig[1] - mu.max)/sd.max,
                                                     -(eig[2] - mu.min)/sd.min)
  return(t)
}

Bootstrap_ERGM = function(hatA, n, numBoot, Phat, MLE.fit){
  eig <- c(eigen(hatA)$values[1], eigen(hatA)$values[n])
  eig.btsp <- matrix(rep(0, numBoot*2), ncol = numBoot)
  
  
  for (i in 1:numBoot) {
    G.btsp <- unname(as.matrix.network.adjacency(simulate(MLE.fit, nsim = 1, burnin = 1000, interval = 1000)))
    A.btsp = (G.btsp - Phat) / sqrt( (n - 1) * Phat * (1 - Phat)) 
    diag(A.btsp) = 0
    
    eig.btsp[1, i] = eigen(A.btsp)$values[1]
    eig.btsp[2, i] = eigen(A.btsp)$values[n]  
  }
  
  mu.max = mean(eig.btsp[1,])
  mu.min = mean(eig.btsp[2,])
  sd.max = sd(eig.btsp[1,])
  sd.min = sd(eig.btsp[2,])
  t = -1.2065335745820 + sqrt(1.607781034581) * max( (eig[1] - mu.max)/sd.max,
                                                     -(eig[2] - mu.min)/sd.min)
  return(t)
}