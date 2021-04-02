# For symmetric Graph G, return its bootstrap statistics
# using MLE. Doing Bootstrap for 50 times. 

return_bootstrap_statistic <- function (n, df, model) {
  
  # Number of bootstraps
  btsp <- 50
  
  if (model != "ER") {
    # Estimated Phat
    Phat <- matrix(unlist(df$Phat), ncol = n)
  }
  
  # Eigenvalues
  eig <- unlist(df$eig)
  
  eig.btsp <- matrix(rep(0, btsp*2), ncol = btsp)
  
  if (model == "SBM" || model == "Beta") {
    for (j in 1:btsp) {
      G.btsp = matrix(rep(0, n^2), ncol = n) 
      G.btsp[upper.tri(G.btsp)] = runif(choose(n, 2)) < Phat[upper.tri(Phat)]
      G.btsp = G.btsp + t(G.btsp)
      
      A.btsp = (G.btsp - Phat) / sqrt( (n - 1) * Phat * (1 - Phat)) 
      diag(A.btsp) = 0
      
      eig.btsp[1, j] = eigen(A.btsp)$values[1]; eig.btsp[2, j] = eigen(A.btsp)$values[n]
    }
  }
  
  if (model == "ER") {
    phat = df$MLE
    for (j in 1:btsp) {
      G.btsp = matrix(rep(0, n^2), ncol = n) 
      G.btsp[upper.tri(G.btsp)] = runif(choose(n, 2)) < phat
      G.btsp = G.btsp + t(G.btsp)
      
      A.btsp = (G.btsp - phat) / sqrt( (n - 1) * phat * (1 - phat)) 
      diag(A.btsp) = 0
      
      eig.btsp[1, j] = eigen(A.btsp)$values[1]; eig.btsp[2, j] = eigen(A.btsp)$values[n]
    }
  }
  
  if (model == "DCSBM") {
    for (i in 1:btsp) {
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
      diag(A.btsp) = rep(0, n)
      
      eig.btsp[1, i] = eigen(A.btsp)$values[1]; eig.btsp[2, i] = eigen(A.btsp)$values[n]
    }
  }
  
  if (model == "ERGM") {
    # MLE of ERGMs
    coef = unlist(df$MLE)
    
    # Formula of MLE
    formula = df$formula
    
    for (j in 1:btsp) {
      g.btsp <- eval(parse(text = paste("simulate(network(n, directed = F) ~", formula, ", coef = coef, directed = F)")))
      G.btsp = unname(as.matrix.network.adjacency(g.btsp))
      
      A.btsp = (G.btsp - Phat) / sqrt( (n - 1) * Phat * (1 - Phat)) 
      diag(A.btsp) = 0
      
      eig.btsp[1, j] = eigen(A.btsp)$values[1]; eig.btsp[2, j] = eigen(A.btsp)$values[n]
    }
  }
  
  if (model == "SRM") {
    symm = df$symm
    if (is.null(symm)) {symm = T}
    
    if (symm == T) {
      ahat = unlist(df$ahat)
      muhat = df$muhat
      ehat = df$ehat
      
      temprhat <- matrix( rep(ahat, n), ncol = n, byrow = F)
      tempchat <- matrix( rep(ahat, n), ncol = n, byrow = T)
      Zhat = muhat + temprhat + tempchat
      
      for (j in 1:btsp) {
        Ehat.vec <- rnorm(choose(n, 2), 0, ehat)
        Ehat <- matrix(rep(0, n^2), ncol = n)
        Ehat[upper.tri(Ehat, diag = F)] = Ehat.vec
        Ehat = Ehat + t(Ehat)
        G.btsp = 1 * (Zhat + Ehat > 0)
        diag(G.btsp) = rep(0, n)
        
        A.btsp = (G.btsp - Phat) / sqrt( (n - 1) * Phat * (1 - Phat)) 
        diag(A.btsp) = 0
        
        eig.btsp[1, j] = eigen(A.btsp)$values[1]; eig.btsp[2, j] = eigen(A.btsp)$values[n]
      }
      
    } else {
      ahat = unlist(df$ahat)
      bhat = unlist(df$bhat)
      muhat = df$muhat
      ehat = df$ehat
      rhohat = df$rhohat

      temprhat <- matrix( rep(ahat, n), ncol = n, byrow = F)
      tempchat <- matrix( rep(bhat, n), ncol = n, byrow = T)
      Zhat <- muhat + temprhat + tempchat

      for (i in 1:btsp) {
        Ehat <- matrix(rep(0, n^2), ncol = n)
        for (j in 1: (n-1)) {
          for (k in (j+1): n) {
            e <- rmvnorm(1, rep(0, 2), matrix(c(1, rhohat, rhohat, 1), ncol = 2, byrow = T) * ehat)
            Ehat[j, k] = e[1]; Ehat[k, j] = e[2]
          }
        }
        G.btsp = 1 * (Zhat + Ehat > 0)
        diag(G.btsp) = rep(0, n)
        
        A.btsp = (G.btsp - Phat) / sqrt( Phat * (1 - Phat))
        diag(A.btsp) = rep(0, n)
        
        eig.btsp[1, i] = eigen(A.btsp %*% t(A.btsp))$values[1]
        eig.btsp[2, i] = eigen(A.btsp %*% t(A.btsp))$values[n]
      }
      
      mu.max = mean(eig.btsp[1,])
      mu.min = mean(eig.btsp[2,])
      sd.max = sd(eig.btsp[1,])
      sd.min = sd(eig.btsp[2,])
      Test = 1 + max( (eig[1] - mu.max)/sd.max, -(eig[2] - mu.min)/sd.min)
      return (Test)
    }
  }
  mu.max = mean(eig.btsp[1,])
  mu.min = mean(eig.btsp[2,])
  sd.max = sd(eig.btsp[1,])
  sd.min = sd(eig.btsp[2,])
  Test = -1.2065335745820 + sqrt(1.607781034581) * max( (eig[1] - mu.max)/sd.max,
                                                        -(eig[2] - mu.min)/sd.min)
  return (Test)
}

