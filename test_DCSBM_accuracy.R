# Type I error test for DCSBM
library(purrr)
library(randnet)

n <- 50

m <- 2

m.vec <- rdunif(n, m)

B.vec <- runif(choose(m, 2), 0, 1)
B <- matrix(rep(0, m^2), ncol = m)
B[upper.tri(B)] = B.vec
B = B + t(B)
diag(B) = runif(m, 0, 1)

theta <- runif(n, 0, 5)

# What if the graph has self loops?
P <- matrix(rep(0, n^2), ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      P[i, j] = B[m.vec[i], m.vec[j]] * theta[i] * theta[j]
    }
  }
}
numSim <- 200
result <- rep(0, 200)
for (index in 1:numSim) {
  print(index)
  Phat <- matrix(rep(0, n^2), ncol = n)
  while (length(which(Phat == 0)) != n) {
    G <- matrix(rep(0, n^2), ncol = n) 
    for (i in 1:n) {
      for (j in 1:n) {
        if (i < j) {
          G[i, j] = rpois(1, P[i, j])
        }
      }
    }
    G = G + t(G)
    
    MLE <- DCSBM.estimate(G, m.vec)
    Phat <- MLE$Phat
    diag(Phat) = rep(0, n)
  }
  Ahat = (G - Phat) / sqrt( (n - 1)* Phat )
  diag(Ahat) = rep(0, n)
  eig <- c(eigen(Ahat)$values[1], eigen(Ahat)$value[n])
  
  # Bootstrap
  B <- 50
  eig.btsp <- matrix(rep(0, B * 2), ncol = B)
  for (i in 1:B) {
    G.btsp <- matrix(rep(0, n^2), ncol = n) 
    for (j in 1:n) {
      for (k in 1:n) {
        if (j < k) {
          G.btsp[j, k] <- rpois(1, Phat[j, k])
        }
      }
    }
    G.btsp = G.btsp + t(G.btsp)
    
    A.btsp = (G.btsp - Phat) / sqrt( (n - 1) * Phat)
    diag(A.btsp) = 0
    
    eig.btsp[1, i] = eigen(A.btsp)$values[1]; eig.btsp[2, i] = eigen(A.btsp)$values[n]
  }
  
  mu.max = mean(eig.btsp[1,])
  mu.min = mean(eig.btsp[2,])
  sd.max = sd(eig.btsp[1,])
  sd.min = sd(eig.btsp[2,])
  result[index] = -1.2065335745820 + sqrt(1.607781034581) * max( (eig[1] - mu.max)/sd.max,
                                                                 -(eig[2] - mu.min)/sd.min)
}

result 
result <- result[1:104]

hist(result, breaks = 30, freq = F)
max(result)
