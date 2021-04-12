# Take a data.frame as input, generate a graph matrix via Beta model, or
# take in a graph matrix, and do a goodness of fit test with respect to 
# Beta model.

generate_Beta <- function (df) {
  
  # Number of nodes
  n <- df$node
  
  # link function
  if (df$link == "exp") {link <- function(temp) {return (exp(-temp))}}
  if (df$link == "expit") {link <- function(temp) {return (exp(temp) / (1 + exp(temp)))}}
  
  # GOF indicator
  fit <- df$fit
  
  if (is.null(fit)) {fit = F}
  
  if (fit == T) {
    # Input G matrix
    G <- matrix(unlist(df$Graph), ncol = n)
    
    # MLE 
    # Now compute MLE with fixed point method. 
    r = function(x, i, j){
      return(1/(exp(-x[j]) + exp(x[i])))
    }
    
    deg = rowSums(G) # compute the degree of G 
    
    phi = function(x){
      x_next = rep(0, n)
      for(i in 1:n){
        temp_sum = 0 
        for(j in 1:n){
          if(i != j){
            temp_sum = temp_sum + r(x, i, j)
          }
        }
        x_next[i] = log(deg[i]) - log(temp_sum)
      }
      return(x_next)
    }
    
    T = 50
    x = runif(n, 0, 1)
    for(index in 1:T){
      x = phi(x) 
    }
    
    betaHat = x 
    
    temphat1 <- matrix(rep(betaHat, n), ncol = n, byrow = T)
    temphat2 <- matrix(rep(betaHat, n), ncol = n, byrow = F)
    temphat <- temphat1 + temphat2
    Phat = link(temphat)
    diag(Phat) = rep(0, n)
    
    Ahat <- (G - Phat) / sqrt( (n - 1) * Phat * (1 - Phat))
    diag(Ahat) = rep(0, n)
    eig <- c(eigen(Ahat)$values[1], eigen(Ahat)$values[n])
    
    result <- data.frame(MLE = I(list(betaHat)), eig = I(list(eig)), Phat = I(list(Phat)))
    return (result)
    
  } else {
    
    # Input Beta parameters, of length n
    beta <- unlist(df$beta)
    
    # Probability matrix P generated via beta model
    temp1 <- matrix(rep(beta, n), ncol = n, byrow = T)
    temp2 <- matrix(rep(beta, n), ncol = n, byrow = F)
    temp <- temp1 + temp2
    
    P <- link(temp)
    
    G = matrix(rep(0, n^2), ncol = n)
  
    
    G[upper.tri(G)] = runif(choose(n, 2)) < P[upper.tri(P)]
    G = G + t(G)
    
    return (G)
  }
}