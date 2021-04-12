alpha <- 0.05

n <- 150

test_model = "Beta"
i = 1 + n - 150

source("Run_Simulations.R")

df_test<- data.frame(
  node = n,
  # Remove the link if test model is not beta.
  link <- "expit",
  fit = T
)
reject <- rep(2, 3)

# true model == ER
file <- read.csv("ER_p.csv")
p <- file[1, i]

df_true = data.frame(
  node = n,
  p <- p
)

reject[1] = run_simulations(n, df_true, df_test, "ER", test_model, alpha)

#true model == SBM

k <- 3

m.file <- read.csv(paste("SBM_m_n", n, ".csv", sep = ""))
m <- unlist(m.file[i, ])

B.file <- read.csv("SBM_B.csv")
B.vec <-unlist(B.file[i, ])
B <- matrix(rep(0, k^2), ncol = k)
B[upper.tri(B)] = B.vec
B = B + t(B)
B.diag.file = read.csv("SBM_B_diag.csv")
B.diag.vec = unlist(B.diag.file[i, ])
diag(B) = B.diag.vec 

df_true = data.frame(
  node = n,
  k = 3,
  m = I(list(m)),
  B = I(list(B))
)

reject[2] = run_simulations(n, df_true, df_test, "SBM", test_model, alpha)


# true model == beta
file <- read.csv(paste("beta_n", n, ".csv", sep = ""))
vec <- t(t(file[i, ]))

df_true = data.frame(
  node = n,
  beta <- I(list(vec)),
  link <- "expit"
)
reject[3] = run_simulations(n, df_true, df_test, "Beta", test_model, alpha)

# write.csv(reject, paste("Data/result_true_", true_model,"_test_", test_model, "_n_", n, "_beta_vec_", node, "_sim_", num, ".csv", sep = ""), row.names = FALSE)
