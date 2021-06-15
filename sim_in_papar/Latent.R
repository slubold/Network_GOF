# LS estimation using projected gradient descent in Ma et al. 2020
# Without covariant matrix X
library(RMTstat)

# Number of iterations 
n.vec = c(100, 200)
# True dimension
d_true.v = c(1, 2, 3, 4)
# Hypothesis dimension
d.fit = c(1, 2, 3, 4)
# Number of sets of parameters
simul = 50
# Final result
final.res <- matrix(rep(0, 2 * simul * 4), nrow = 2)
for (out in 1:2) {
  n = n.vec[out]
  J = diag(n) - (1 / n) * rep(1, n) %*% t(rep(1, n))
  for (inn in 1:4) {
    d_true = d_true.v[inn]
    
    for (numSim in 1:simul) {
      # True values
      Z_true <- mvtnorm::rmvnorm(n, mean = rep(0, d_true), sigma = diag(d_true))
      a_true = runif(n, -0.02, -0.01)
      theta_true <- a_true %*% t(rep(1, n)) + rep(1, n) %*% t(a_true) + Z_true %*% t(Z_true)
      
      # Sigmoid function  
      sigmoid <- function(x) {exp(x)/(1+exp(x))} 
      
      # P matrix
      P = sigmoid(theta_true); diag(P) = rep(0, n)
      tw.stat <- rep(0, 100)
      for (num in 1:100) {
        result <- 0
        inner = 0
        while ((inner <= 3) && (result == 0)) {
          inner = inner + 1
          Ahat <- matrix(rep(0, n^2), ncol = n)
          while ( (length(which(Ahat == 0)) > n) || (length(which(Ahat == 1)) > 0) ) {
            A = matrix(rep(0, n^2), ncol = n); A[upper.tri(A)] = runif(choose(n, 2)) < P[upper.tri(P)]; A = A + t(A)
            # Initializing: SVT (Robust)
            p_ave = sum(A)/(n*(n-1))
            tau0 = sqrt(n * p_ave)
            temp = svd(A)
            d = diag(temp$d) %*% (diag(temp$d) > tau0)
            temp = temp$u %*% d %*% t(temp$v)
            M = 4 # By default
            temp = pmax(temp, 1/(1 + exp(M)))
            temp = pmin(temp, 1/(1 + exp(-M)))
            logit = -log(1 / temp - 1)
            logit = (logit + t(logit)) / 2
            logitsum = rowSums(logit)
            alpha0 = (logitsum - 1 / 2  * rep(1, n) * mean(logitsum)) / n;
            R = logit - alpha0 %*% t(rep(1, n)) - rep(1, n) %*% t(alpha0)
            R = R - rep(1, n) %*% (rep(1, n) %*% R / n)
            G0 = R - (R / n) %*% rep(1, n) %*% t(rep(1, n))
            G0 = (G0 + t(G0)) / 2
            eig = eigen(G0)
            d = d.fit[inner]
            value <- sort(eig$values, decreasing = T)[1:d]
            indice <- c()
            for (i in 1:d) {  
              indice <- c(indice, which(eig$values == value[i]))
            }
            if (d != 1) {
              Z = eig$vectors[, indice] %*% sqrt(diag(value))
            } else {
              Z = eig$vectors[, indice] * sqrt(value)
            }
            a = alpha0
            eta = 2
            zstep = 1 / ((norm(Z, type = "2"))^2)
            astep = 2 / n#eta / (2 * n)
            # Theta matrix
            theta <- a %*% t(rep(1, n)) + rep(1, n) %*% t(a) + Z %*% t(Z)
            change_alpha = 1.0
            change_G = 1.0
            maxiter = 500
            miniter = 100
            # Thershold
            epsilon = 1e-5
            
            iter = 0; report_every = 100
            while (((max(abs(change_alpha), abs(change_G)) > epsilon) && (iter<maxiter))||(iter<miniter)) {
              iter = iter + 1
              temp.Z = Z + 2 * zstep * (A - sigmoid(theta)) %*% Z
              temp.Z = J %*% temp.Z
              temp.a = a + 2 * astep * (A - sigmoid(theta)) %*% rep(1, n)
              if (iter %% report_every == 0) {
                G0 = Z %*% t(Z); G1 = temp.Z %*% t(temp.Z)
                change_alpha = norm(temp.a - a, "2") / norm(a, "2")
                change_G = norm(G1 - G0, "f")/norm(G0, "f")
              }
              # Update
              Z = temp.Z; a = temp.a; theta <- a %*% t(rep(1, n)) + rep(1, n) %*% t(a) + Z %*% t(Z)
            }
            
            Phat <- sigmoid(a %*% t(rep(1, n)) + rep(1, n) %*% t(a) + Z %*% t(Z)); diag(P) = rep(0, n)
            Ahat = (A - Phat) / sqrt((n - 1) * Phat * (1 - Phat)); diag(Ahat) = rep(0, n)
          }
          eig.observed = eigen(Ahat)$values
          # Bootstrap
          B = 50
          eig.btsp = matrix(rep(0, B * 2), ncol = 2)
          for(indexBoot in 1:B){
            G.sim = matrix(rep(0, n^2), ncol = n)
            G.sim[upper.tri(G.sim)] = runif(choose(n, 2)) < Phat[upper.tri(Phat)]
            G.sim = G.sim + t(G.sim)
            A.sim = (G.sim - Phat)/sqrt((n - 1) * Phat * (1- Phat))
            diag(A.sim) = rep(0, n)
            eig.btsp[indexBoot,] = c(eigen(A.sim)$values[1], eigen(A.sim)$values[n])  
          }
          
          mu.max = mean(eig.btsp[,1])
          mu.min = mean(eig.btsp[,2])
          sd.max = sd(eig.btsp[,1])
          sd.min = sd(eig.btsp[,2])
          test = -1.2065335745820 + sqrt(1.607781034581) * max( (eig.observed[1] - mu.max)/sd.max,
                                                                -(eig.observed[n] - mu.min)/sd.min)
          result <- (test < qtw(1 - (0.05 / 2))) && (test > qtw(0.05 / 2)) * 1
        }
        if (result == 1) {
          if (inner == d_true) {
            tw.stat[num] = 1
          }
        }
      }
      final.res[out, (inn - 1) * simul + numSim] = mean(tw.stat)
    }
  }
}
print(final.res)
