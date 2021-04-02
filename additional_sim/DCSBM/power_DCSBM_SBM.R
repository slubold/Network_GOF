# Suppose community number k = 3 (Remember to change this also in the generating file!)

library(randnet)

k = 3

alpha <- 0.05

# Manually adjust n
n = 150
temp = n - 150

true_model = "DSBM"
test_model = "DCSBM"

source("Run_Simulations.R")

# G <- generate_DSBM(df_true)
# result <- generate_DSBM(df_test)

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


df_test<- data.frame(
  node = n,
  k = 3,
  m = I(list(m)),
  fit = T
)

reject = run_simulations(n, df_true, df_test, true_model, test_model, alpha)

# We just want one simulation: 
write.csv(reject, paste("Data/result_true_", true_model,"_test_", test_model, "_n_", n, 
                        "_DCSBM_", n, "_sim_", num, ".csv", sep = ""), row.names = FALSE)

