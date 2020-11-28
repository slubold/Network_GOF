# Take a data.frame as input, generate a graph matrix via SBM model, or
# take in a graph matrix, and do a goodness of fit test with respect to 
# SBM model.

generate_SBM <- function (df) {
  
  # Number of nodes
  n <- df$node

  # GOF indicator
  fit <- df$fit
  
  if (is.null(fit)) {fit = F}
  
  if (fit == T) {
    
    # Input G matrix
    G <- matrix(unlist(df$Graph), ncol = n)
    g = graph_from_adjacency_matrix(G, mode = "undirected", diag = F)
    
    # greedy method (hiearchical, fast method)
    c1 = cluster_leading_eigen(g) # = cluster_fast_greedy(g)
    
    # memberships of nodes
    m <- membership(c1)
    
    # number of communitys m
    k <- length(c1)
    
    # sum of total nodes in each community
    sum.vec <- rep(0, k)
    for (i in 1:n) {
      sum.vec[m[i]] = sum.vec[m[i]] + 1
    }
    
    # Coefficient of Bhat
    temp1 <- matrix(rep(sum.vec, k), ncol = k, byrow = T)
    temp2 <- matrix(rep(sum.vec, k), ncol = k, byrow = F)
    coef <- temp1 * temp2
    diag(coef) = diag(coef) * (sum.vec - 1)/ 2 /sum.vec
    
    # Compute MLE with respect to G
    # Bhat for different communities
    sum.node <- matrix(rep(0, k^2), ncol = k)
    for (j in 1:n) {
      for (l in 1:n) {
        if (G[j, l] == 1) {
          sum.node[m[j], m[l]] = sum.node[m[j], m[l]] + 1
        }
      }
    }
    diag(sum.node) = diag(sum.node) / 2
    
    Bhat <- sum.node / coef
    
    Phat <- matrix(rep(0, n^2), ncol = n)
    for (l in 1:n) {
      for (j in 1:n) {
        Phat[l, j] = Bhat[m[l], m[j]]
      }
    }
    diag(Phat) = rep(0, n)
    
    Ahat = (G - Phat) / sqrt( (n - 1) * Phat * (1 - Phat)) 
    diag(Ahat) = rep(0, n)
    
    eig <- c(eigen(Ahat)$values[1], eigen(Ahat)$values[n])
    
    result <- data.frame(MLE = I(list(Bhat)), eig = I(list(eig)), Phat = I(list(Phat)))
    
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