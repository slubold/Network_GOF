library(RMTstat)

# Parameters in ARD sections
gamma = 1/3
n = 150
beta = 1
m = gamma * n
numSim = 200
test = rep(0, numSim)

# ARD for ER model
return_value = function(n, null, numSim){
  
  pstar = runif(1, 0.2, 0.8) # parameter for ER model
  gamma_m = 1/2
  gamma_k = 1/10
  m = gamma_m * n
  K = gamma_k * n
  membership = rep(1:K, each=1/gamma_k)
  test = rep(0, n)
  test.complete = rep(0, n)
  beta = 1
  
  for(index in 1:numSim){
    
    # Generate ER model using pstar
    if(null ==  TRUE){
      g = igraph::erdos.renyi.game(n, pstar)
    }
    
    # Two community SBM, equal number
    if(null == FALSE){
      g = igraph::sbm.game(n, matrix(c(1/10, 1/300, 1/300, 1/10), ncol = 2), block.sizes =  c(n/2, n/2))
    }
    
    Adj = matrix(igraph::get.adjacency(g), ncol = n)
    
    # MLE of p, given true matrix
    pHat.complete = sum(Adj)/(n*(n - 1)) 
    
    A.complete = (Adj - pHat.complete) / sqrt((n - 1) * pHat.complete * (1 - pHat.complete))
    diag(A.complete) = rep(0, n)
    
    # Tracy-Widom statistic, Thm 1
    test.complete[index] = n^(2/3) * (eigen(A.complete)$values[1] - 2)
    
    # ARD 
    Y = matrix(rep(0, m * K), ncol = K)
    for(i in 1:m){
      for(j in 1:K){
        Y[i,j] = sum(g[i, membership == j])
      }
    }
    # MLE of ARD
    pHat = mean(Y/10)
    
    # A = (Y - pstar * 1/gamma_k)/sqrt(pstar * 1/gamma_k * (1 - pstar))
    A = (Y - pHat * (1 / gamma_k)) / sqrt(pHat * 1/gamma_k * (1 - pHat))

    mu = (sqrt(m + beta - 2) + sqrt(K))^2
    sigma = sqrt(mu) * (1/sqrt(m + beta - 2) + 1/sqrt(K))^(1/3)
    test[index] = (eigen(t(A) %*% A)$values[1] - mu)/sigma
  }
  rej.ARD = mean(test > RMTstat::qtw(0.975, beta)) + mean(test < RMTstat::qtw(0.025, beta))
  rej.complete = mean(test.complete > RMTstat::qtw(0.975)) + mean(test.complete < RMTstat::qtw(0.025))
  return(c(rej.ARD, rej.complete))
}


n.vec = c(30, 60, 90)
numOuter = 25

values.true.ARD = matrix(rep(0, length(n.vec) * numOuter), ncol = numOuter)
values.true.complete = matrix(rep(0, length(n.vec) * numOuter), ncol = numOuter)

for(i in 1:length(n.vec)){
  for(j in 1:numOuter){
    temp = return_value(n.vec[i], TRUE, 100)
    values.true.complete[i,j] = temp[2]
    values.true.ARD[i,j] = temp[1]
  }
}

values.F.ARD = matrix(rep(0, length(n.vec) * numOuter), ncol = numOuter)
values.F.complete = matrix(rep(0, length(n.vec) * numOuter), ncol = numOuter)

for(i in 1:length(n.vec)){
  for(j in 1:numOuter){
    temp = return_value(n.vec[i], FALSE, 100)
    values.F.complete[i,j] = temp[2]
    values.F.ARD[i,j] = temp[1]
  }
}

write.csv(values.true.complete, "er_error_complete.csv")
write.csv(values.true.ARD, "er_error_ARD.csv")
write.csv(values.F.complete, "sbm_er_power_complete.csv")
write.csv(values.F.ARD, "sbm_er_power_ARD.csv")

