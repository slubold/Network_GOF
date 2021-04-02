directed_func <- function(n, null, numSim){
  
  # ER model
  if (null == "true"){
    p = runif(1)
    P = matrix(rep(p, n^2), ncol = n)
  }
  
  # SBM model
  if (null == "false"){
    P = matrix(rep(0, n^2), ncol = n)
    P.group = matrix(rep(0, 4), ncol = 2)
    mem = c(rep(1, n/2), rep(2, n/2))
    P.group[1,1] = runif(1, 0, 0.2)
    P.group[1,2] = runif(1, 0, 0.2)
    P.group[2,1] = runif(1, 0, 0.2)
    P.group[2,2] = runif(1, 0, 0.2)
    P[mem == 1, mem == 1] = P.group[1,1]
    P[mem == 1, mem == 2] = P.group[1,2]
    P[mem == 2, mem == 1] = P.group[2,1]
    P[mem == 2, mem == 2] = P.group[2,2]
  }
  
  test = rep(0, numSim)
  test.TW2 = rep(0, numSim)
  test.boot = rep(0, numSim)
  
  for(i in 1:numSim){
    G = matrix(as.numeric(runif(n^2, 0, 1)) < P, ncol = n)
    diag(G) = rep(0, n)
    
    pHat = sum(G)/(n^2)
    
    A = (G - pHat)/sqrt(pHat * (1 - pHat))
    diag(A) = rep(0, n)
    
    test[i] = eigen(A %*% t(A))$values[n] * n  
    test.TW2[i] = (eigen(A %*% t(A))$values[1] - 4 * n)/(2 * sqrt(n)  * (2 * sqrt(n))^(-1/3))
    
    # Do bootstrapping 
    M = 100
    lambda.boot = matrix(rep(0, M ), ncol = 1)
    for(index.boot in 1:M){
      A.boot = matrix(rep(0, n^2), ncol = n)
      A.boot=  matrix(as.numeric(runif(n^2, 0, 1)) < pHat, ncol = n)
      B.boot = (A.boot - pHat)/sqrt( (1 - pHat) * pHat)
      diag(B.boot) = rep(0, n)
      lambda.boot[index.boot] = (eigen(B.boot %*% t(B.boot))$values[1])
    }
    mu.1 = mean(lambda.boot)
    s.1 = sd(lambda.boot)
    mu.tw = -1.2065335745820
    s.tw = sqrt(1.607781034581)
    test.boot[i] = mu.tw + s.tw * (eigen(A %*% t(A))$values[1] - mu.1)/s.1
  }
  return(c(test.boot, test, test.TW2))
}


numOuter = 25
numSim = 100
null = "false"
n_vec = c(25, 50, 100)
reject.boot = matrix(rep(0, length(n_vec) * numOuter), ncol = numOuter)
reject.exp = matrix(rep(0, length(n_vec) * numOuter), ncol = numOuter)
reject.TW = matrix(rep(0, length(n_vec) * numOuter), ncol = numOuter)

for(indexN in 1:length(n_vec)){
  for(indexOuter in 1:numOuter){
    temp = directed_func(n_vec[indexN], null, numSim)
    temp.boot = temp[1:numSim]
    temp.exp = temp[seq(from = numSim +1, to = 2 * numSim, by = 1)]
    
    temp.TW = temp[seq(from = 2*numSim +1, to = 3 * numSim, by = 1)]
    
    
    reject.boot[indexN, indexOuter] = mean(temp.boot > RMTstat::qtw(0.975)) + mean(temp.boot < RMTstat::qtw(0.025))
    reject.TW[indexN, indexOuter] = mean(temp.TW > RMTstat::qtw(0.975)) + mean(temp.TW < RMTstat::qtw(0.025))
    reject.exp[indexN, indexOuter] = mean(temp.exp > qexp(0.95))
  }
}