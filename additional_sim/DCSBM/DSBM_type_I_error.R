# DCSBM_weighted type I error test
# Suppose community number k = 3 (Remember to change this also in the generating file!)

library(randnet)
library(RMTstat)

k = 3

alpha <- 0.05

# Manually adjust n
n = 150
temp = n - 150

true_model = "DSBM"
test_model = "DSBM"

source("Generate_DSBM.R")
source("Return_Bootstrap_Statistic.R")
# Suppose args start at 1
# args = commandArgs(TRUE)
# Manually adjust i

i = 1
m.file <- read.csv(paste("DCSBM_m_n", n, ".csv", sep = ""))
m <- unlist(m.file[i, ])

B.file <- read.csv("DCSBM_B.csv")
B.vec <-unlist(B.file[i + temp, ])
B <- matrix(rep(0, k^2), ncol = k)
B[upper.tri(B)] = B.vec
B = B + t(B)
B.diag.file = read.csv("DCSBM_B_diag.csv")
B.diag.vec = unlist(B.diag.file[i + temp, ])
diag(B) = B.diag.vec 

df_true = data.frame(
  node = n,
  k = 3,
  m = I(list(m)),
  B = I(list(B))
)

Phat <- matrix(rep(0, n^2), ncol = n)
while (length(which(Phat == 0)) != n) {
  G <- generate_DSBM(df_true)
  Phat <- SBM.estimate(G, m)$Phat
  diag(Phat) = rep(0, n)
}

df_test<- data.frame(
  Graphs <- I(list(G)),
  node = n,
  k = 3,
  m = I(list(m)),
  fit = T
)

result <- generate_DSBM(df_test)
test <- return_bootstrap_statistic(n, result, test_model)
reject = 1
if ( (test < qtw(1 - (alpha / 2))) && (test > qtw(alpha / 2)) ) {reject = 0}

# We just want one simulation: 
write.csv(reject, paste("Data/result_true_", true_model,"_test_", test_model, "_n_", n, 
                        "_DSBM_", n, "_sim_", num, ".csv", sep = ""), row.names = FALSE)

