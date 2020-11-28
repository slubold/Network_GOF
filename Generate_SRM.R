# Take a data.frame as input, generate a graph matrix via SRM model, or
# take in a graph matrix, and do a goodness of fit test with respect to 
# SRM model.

generate_SRM <- function (df) {
  
  # Number of nodes
  n <- df$node
  
  symm <- df$symm
  if (is.null(symm)) {symm = T}
  
  # GOF indicator
  fit <- df$fit
  
  if (is.null(fit)) {fit = F}
  
  if (fit == T) {
    
    Y = matrix(unlist(df$Graph), ncol = n)
    if (symm == F) {
      
      MLE.fit <- ame(Y, model = "bin")
      ahat <- MLE.fit[["APM"]]
      bhat <- MLE.fit[["BPM"]]
      temprhat <- matrix( rep(ahat, n), ncol = n, byrow = F)
      tempchat <- matrix( rep(bhat, n), ncol = n, byrow = T)
      muhat <- mean(MLE.fit$BETA) 
      Zhat <- muhat + temprhat + tempchat
      rhohat <- mean(MLE.fit$VC[,4])
      ehat <- mean(MLE.fit$VC[,5])
      
      numSim <- 200
      Phat <- matrix(rep(0, n^2), ncol = n)
      for (i in 1:numSim) {
        Ehat <- matrix(rep(0, n^2), ncol = n)
        for (j in 1: (n-1)) {
          for (k in (j+1): n) {
            e <- rmvnorm(1, rep(0, 2), matrix(c(1, rhohat, rhohat, 1), ncol = 2, byrow = T) * ehat)
            Ehat[j, k] = e[1]; Ehat[k, j] = e[2]
          }
        }
        Yhat = 1 * (Zhat + Ehat > 0)
        diag(Yhat) = rep(0, n)
        Phat = Phat + Yhat
      }
      Phat <- Phat / numSim
      Phat.vec = pmax(0.01, Phat)
      Phat.vec = pmin(0.99, Phat.vec)
      Phat <- matrix(Phat.vec, ncol = n, byrow = F)
      diag(Phat)= rep(0, n)
      
      A = (Y - Phat) / sqrt( Phat * (1 - Phat))
      diag(A) = rep(0, n)
      
      test <- n * eigen(A %*% t(A))$values[n]
      result <- data.frame(test = test, 
                           ahat = I(list(ahat)),
                           bhat = I(list(bhat)),
                           muhat = muhat,
                           rhohat = rhohat,
                           ehat = ehat,
                           Phat = I(list(Phat)),
                           symm = F
      )
      return (result)
    }
    
    MLE.fit <- ame(Y, model = "bin", symmetric = TRUE)
    ahat <- MLE.fit[["APM"]]
    temprhat <- matrix( rep(ahat, n), ncol = n, byrow = F)
    tempchat <- matrix( rep(ahat, n), ncol = n, byrow = T)
    muhat <- mean(MLE.fit$BETA) 
    ehat <- mean(MLE.fit$VC[,2])
    Zhat = muhat + temprhat + tempchat
    
    numSim <- 200
    Phat <- matrix(rep(0, n^2), ncol = n)
    for (i in 1:numSim) {
      Ehat.vec <- rnorm(choose(n, 2), 0, ehat)
      Ehat <- matrix(rep(0, n^2), ncol = n)
      Ehat[upper.tri(Ehat, diag = F)] = Ehat.vec
      Ehat = Ehat + t(Ehat)
      Yhat = 1 * (Zhat + Ehat > 0)
      diag(Yhat) = rep(0, n)
      Phat = Phat + Yhat
    }
    Phat <- Phat / numSim
    Phat.vec = pmax(0.01, Phat)
    Phat.vec = pmin(0.99, Phat.vec)
    Phat <- matrix(Phat.vec, ncol = n, byrow = F)
    diag(Phat)= rep(0, n)
    
    A = (Y - Phat) / sqrt( (n - 1) * Phat * (1 - Phat))
    diag(A) = rep(0, n)
    eig <- c(eigen(A)$values[1], eigen(A)$values[n])
    
    result <- data.frame(eig = I(list(eig)),
                         ahat = I(list(ahat)),
                         muhat = muhat,
                         ehat = ehat,
                         Phat = I(list(Phat))
    )
    
    return (result)
    
  } else {
    
    mu <- df$mu
    a.vec <- unlist(df$a.vec)
    sigma.e <- df$sigma.e
    
    if (symm == F) {
      rho <- df$rho
      b.vec <- unlist(df$b.vec)
      cov.e <- matrix(c(1, rho, rho, 1), ncol = 2, byrow = T) * sigma.e
      
      Z <- matrix(rep(mu, n^2), ncol = n)
      temprow <- matrix( rep(a.vec, n), ncol = n, byrow = F)
      tempcol <- matrix( rep(b.vec, n), ncol = n, byrow = T)
      Z <- Z + temprow + tempcol 
      E <- matrix(rep(0, n^2), ncol = n)
      for (j in 1: (n-1)) {
        for (k in (j+1): n) {
          e <- rmvnorm(1, rep(0, 2), cov.e)
          E[j, k] = e[1]; E[k, j] = e[2]
        }
      }
      Y = 1 * (Z + E > 0)
      diag(Y) = rep(0, n)
      return(Y)
    }
    
    Z <- matrix(rep(mu, n^2), ncol = n)
    temprow <- matrix( rep(a.vec, n), ncol = n, byrow = F)
    tempcol <- matrix( rep(a.vec, n), ncol = n, byrow = T)
    Z = Z + temprow + tempcol
    
    E.vec <- rnorm(choose(n, 2), 0, sigma.e)
    E <- matrix(rep(0, n^2), ncol = n)
    E[upper.tri(E, diag = F)] = E.vec
    E = E + t(E)
    Y = 1 * (Z + E > 0)
    diag(Y) = rep(0, n)
    return (Y)
  }
}