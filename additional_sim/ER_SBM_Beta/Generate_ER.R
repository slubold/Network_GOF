# Take a data.frame as input, generate a graph matrix via ER model, or
# take in a graph matrix, and do a goodness of fit test with respect to 
# ER model.

generate_ER <- function (df) {
  
  # Number of nodes
  n <- df$node
  
  # GOF indicator
  fit <- df$fit
  
  if (is.null(fit)) {fit = F}
  
  if (fit == T) {
    # Input G matrix
    G <- matrix(unlist(df$Graph), ncol = n)
    
    # MLE 
    phat <- (sum(G)/ 2) / choose(n, 2)
    Ahat <- (G - phat) / sqrt( (n - 1) * phat * (1 - phat))
    diag(Ahat) = rep(0, n)
    eig <- c(eigen(Ahat)$values[1], eigen(Ahat)$values[n])
    
    result <- data.frame(MLE = phat, eig = I(list(eig)), phat = phat)
    return (result)
  } else {
    
    # probability of ER model
    p <- df$p
    
    # Draw the random graph G
    G <- matrix(rep(0, n^2), ncol = n) 
    G[upper.tri(G)] <- runif(choose(n, 2)) < p  
    G = G + t(G)
    return (G)
  }
}