generate_DSBM <- function (df) {
  
  # Number of nodes
  n <- df$node
  
  # GOF indicator
  fit <- df$fit
  
  if (is.null(fit)) {fit = F}
  
  if (fit == T) {
    
    # Input G matrix
    G <- matrix(unlist(df$Graph), ncol = n)
    # g = graph_from_adjacency_matrix(G, mode = "undirected", diag = F)
    
    # greedy method (hiearchical, fast method)
    # c1 = cluster_leading_eigen(g) # = cluster_fast_greedy(g)
    
    # memberships of nodes
    m <- unlist(df$m)
    
    # number of communitys m
    k <- df$k
    
    MLE <- SBM.estimate(G, m)
    Phat <- MLE$Phat
    diag(Phat) = rep(0, n)
    
    Ahat = (G - Phat) / sqrt( (n - 1) * Phat) 
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
      if (i != j) {
        P[i, j] = B[m[i], m[j]]
      }
    }
  }
  
  Phat <- matrix(rep(0, n^2), ncol = n)
  # while (length(which(Phat == 0)) != n) {
    G <- matrix(rep(0, n^2), ncol = n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i < j) {
          G[i, j] = rpois(1, P[i, j])
        }
      }
    }
    G = G + t(G)
  #   MLE <- DCSBM.estimate(G, m)
  #   Phat <- MLE$Phat
  #   diag(Phat) = rep(0, n)
  # }
  return (G)
}
}