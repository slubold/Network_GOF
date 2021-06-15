# Power test for Beta model with different link function.
# MLE for Beta Model in diaconis paper, using fixed point method
# n: network size
# numSim: Number of simulations to run
# generate: Input link function, "beta_exp", "beta_expit", and "non_par"

library(RMTstat)

# Type 1 error
return_rejection = function(n, numSim, generate, alpha){
  
  p.triangles = rep(0, numSim)
  p.deg = rep(0, numSim)
  Test = rep(0, numSim)
  
  # Generate beta_1, ..., beta_n from Uniform(-2, 0)
  beta.vec <- function (n) {runif(n, -2, 0)}
  beta <- beta.vec(n)
  P = matrix(rep(0, n^2), ncol = n)
  
  # Probability matrix P generated via beta model
  if(generate == "beta_exp"){
    temp1 <- matrix(rep(beta, n), ncol = n, byrow = T)
    temp2 <- matrix(rep(beta, n), ncol = n, byrow = F)
    temp <- temp1 + temp2
    P = exp(temp)
  }
  
  if(generate == "beta_expit"){
    temp1 <- matrix(rep(beta, n), ncol = n, byrow = T)
    temp2 <- matrix(rep(beta, n), ncol = n, byrow = F)
    temp <- temp1 + temp2
    P = exp(temp)/(1+exp(temp))
  }
  
  # Non-parameteric model, correspond to beta parameterization
  if(generate == "non_par"){
    P = matrix(rep(0, n^2), ncol = n)
    P[upper.tri(P)] = runif(choose(n, 2), 0, 0.17)
    P = P + t(P)
    diag(P) = rep(0, n)
  }
  
  Test <- rep(0, numSim)
  
  avg.deg.obs = rep(0, numSim)
  num.tri.obs = rep(0, numSim)
  avg.deg.boot = rep(0, numSim)
  avg.num.tri.boot = rep(0, numSim)
  
  for (indexSim in 1:numSim) {

    # Generate graph G
    G = matrix(rep(0, n^2), ncol = n)
    while(max(rowSums(G) == 0)){
      G = matrix(rep(0, n^2), ncol = n)
      G[upper.tri(G)] = runif(choose(n, 2)) < P[upper.tri(P)]
      G = G + t(G)
    }
    
    avg.deg.obs[indexSim] = mean(rowSums(G))
    num.tri.obs[indexSim] = sum(diag(G %*% G %*% G))/6
    
    # Now compute MLE with fixed point method.
    r = function(x, i,j){
      return(1/(exp(-x[j]) + exp(x[i])))
    }
    
    deg = rowSums(G) # compute the degree of G
    
    phi = function(x){
      x_next = rep(0, n)
      for(i in 1:n){
        temp_sum = 0 # initialize the sum in the definition of phi
        for(j in 1:n){
          if(i != j){
            temp_sum = temp_sum + r(x, i, j)
          }
        }
        x_next[i] = log(deg[i]) - log(temp_sum)
      }
      return(x_next)
    }
    
    T = 30
    x = runif(n, 0, 1)
    for(index in 1:T){
      x = phi(x)
    }
    
    betaHat = x # use the last iteration of the fixed point method as betaHat.
    
    temphat1 <- matrix(rep(betaHat, n), ncol = n, byrow = T)
    temphat2 <- matrix(rep(betaHat, n), ncol = n, byrow = F)
    temphat <- temphat1 + temphat2
    Phat = exp(temphat) / (1 + exp(temphat))
    diag(Phat) = rep(0, n)
    
    Ahat = (G - Phat) / sqrt((n - 1) * Phat * (1 - Phat)) 
    diag(Ahat) = rep(0, n)
    eig.observed = eigen(Ahat)$values
    
    # Now bootstrap
    B = 50
    temp.degree = rep(0, B)
    temp.triangles = rep(0, B)
    eig.btsp = matrix(rep(0, B * 2), ncol = 2)
    for(indexBoot in 1:B){
      G.sim = matrix(rep(0, n^2), ncol = n)
      G.sim[upper.tri(G.sim)] = runif(choose(n, 2)) < Phat[upper.tri(Phat)]
      G.sim = G.sim + t(G.sim)
      A.sim = (G.sim  - Phat)/sqrt((n - 1) * Phat * (1- Phat))
      diag(A.sim) = rep(0, n)
      eig.btsp[indexBoot,] = c(eigen(A.sim)$values[1], eigen(A.sim)$values[n])  
      temp.degree[indexBoot] = mean(rowSums(G.sim))
      temp.triangles[indexBoot] = sum(diag(G.sim %*% G.sim %*% G.sim))/6
    }
    
    mu.max = mean(eig.btsp[,1])
    mu.min = mean(eig.btsp[,2])
    sd.max = sd(eig.btsp[,1])
    sd.min = sd(eig.btsp[,2])
    Test[indexSim] = -1.2065335745820 + sqrt(1.607781034581) * max( (eig.observed[1] - mu.max)/sd.max,
                                                                    -(eig.observed[n] - mu.min)/sd.min)
    
    avg.deg.boot[indexSim] = mean(temp.degree)
    avg.num.tri.boot[indexSim] = mean(temp.triangles)
    
  }
  return(c(Test, avg.deg.obs, num.tri.obs, avg.deg.boot, avg.num.tri.boot))
}


