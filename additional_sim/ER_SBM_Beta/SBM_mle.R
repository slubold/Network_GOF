# SBM MLE 
source("Generate_SBM.R")
source("Generate_ER.R")
source("Generate_Beta.R")

SBM.mle <- function (df) {
  
  model <- df$model
  
  n <- df$node
  
  if (model == "ER" || model == "Beta") {
    
    Phat <- matrix(rep(0, n^2), ncol = n)
  
    while (length(which(Phat == 0)) != n) {
      
      # Input G matrix
      G <- eval(parse(text = paste("generate_", model, "(df)", sep = "")))
      
      g = graph_from_adjacency_matrix(G, mode = "undirected", diag = F)
      
      # greedy method (hiearchical, fast method)
      c1 = cluster_leading_eigen(g)
      
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
    }
  }
  
  if (model == "SBM") {
    m <- unlist(df$m); k <- df$k; B <- matrix(unlist(df$B), ncol = k)
    
    Phat <- matrix(rep(0, n^2), ncol = n)
    
    while (length(which(Phat == 0)) != n) {
      
      # Input G matrix
      G <- generate_SBM(df)
      
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
    }
  }
  
  inform <- data.frame(G= I(list(G)), m = I(list(m)), Phat = I(list(Phat)), k = k)
  return(inform)
}