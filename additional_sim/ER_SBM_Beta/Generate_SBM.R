# Take a data.frame as input, generate a graph matrix via SBM model, or
# take in a graph matrix, and do a goodness of fit test with respect to 
# SBM model.
library(igraph)

generate_SBM <- function (df) {
  
  # Number of nodes
  n <- df$node

  # GOF indicator
  fit <- df$fit
  
  if (is.null(fit)) {fit = F}
  
  if (fit == T) {
    
    Phat <- matrix(unlist(df$Phat), ncol = n)
    
    Ahat = (G - Phat) / sqrt( (n - 1) * Phat * (1 - Phat)) 
    diag(Ahat) = rep(0, n)
    
    eig <- c(eigen(Ahat)$values[1], eigen(Ahat)$values[n])
    
    result <- data.frame(eig = I(list(eig)), Phat = I(list(Phat)))
    
    return (result)
    
  } else {
    
    # Number of community
    k <- df$k
    
    # Community matrix B
    B <- matrix(unlist(df$B), ncol = k)
    
    # Membership vector
    m <- unlist(df$m)
    
    # Probability matrix P
    P <- matrix(rep(0, n^2), ncol = n)
    for (i in 1:n) {
      for (j in 1:n) {
        P[i, j] = B[m[i], m[j]]
      }
    }
    diag(P) = rep(0, n)
    
    # Draw the random graph G
    G <- matrix(rep(0, n^2), ncol = n) 
    G[upper.tri(G)] <- runif(choose(n, 2)) < P[upper.tri(P)]
    G = G + t(G)
    return (G)
  }
}