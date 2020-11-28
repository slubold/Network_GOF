# ERGM model. Generate a ERGM graph by default. Compute the MLE 
# and compute the corresponding eigenvalue otherwise.
library(ergm)
library(sna)
library(RMTstat)

generate_ERGM <- function (df) {
  
  # Number of nodes
  n <- df$node
  
  # Number of simulations
  numSim <- 250
  
  # Formula for ERGM
  formula = df$formula
  
  # GOF indicator
  fit <- df$fit
  
  if (is.null(fit)) {fit = F}
  
  if (fit == T) {
    
    G = matrix(unlist(df$Graph), ncol = n)
    
    g = network(G, directed=FALSE)
    MLE.fit <- eval(parse(text = paste("ergm(g ~ ", formula, ")")))
    
    Phat <- matrix(rep(0, n^2), ncol = n)
    for (i in 1:numSim) {
      # Ergms model 
      ghat = simulate(MLE.fit, nsim = 1, burnin = 1000, interval = 1000) 
      Ghat = unname(as.matrix.network.adjacency(ghat))
      Phat= Phat + Ghat
    }
    Phat <- Phat / numSim # Estimated Phat
    
    Ahat = (G - Phat) / sqrt((n - 1) * Phat * (1 - Phat))
    diag(Ahat) = rep(0, n)
    eig <- c(eigen(Ahat)$values[1], eigen(Ahat)$values[n])
    
    result <- data.frame(MLE = I(list(MLE.fit$coef)), formula = formula,  eig = I(list(eig)), Phat = I(list(Phat)))
    return (result)
    
  } else {
    
    # function of ERGM
    param = unlist(df$param)
    #g = simulate(network(n, directed = F) ~ formula, coef = param, directed = F)
    g <- eval(parse(text = paste("simulate(network(n, directed = F) ~", formula, ",coef = param, directed = F)")))
    G = unname(as.matrix.network.adjacency(g))
    return (G)
  }
  
  
  
  
  
  
  
  
  
}