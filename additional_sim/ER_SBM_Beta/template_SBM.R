library(randnet)
library(RMTstat)

k = 3

alpha <- 0.05

# Manually adjust n
n = 200

source("Run_Simulations.R")
source("SBM_mle.R")
source("Generate_SBM.R")
source("Return_Bootstrap_Statistic.R")

# Suppose args start at 1
# args = commandArgs(TRUE)
# Manually adjust i

i = 1
reject <- rep(1, 3)

# ER
file <- read.csv("ER_p.csv")
p <- file[1, i]

df_true = data.frame(
  node = n,
  p <- p,
  model = "ER"
)

temp <- SBM.mle(df_true)
G <- matrix(unlist(temp$G), ncol = n)
m1 <- unlist(temp$m)
Phat <- matrix(unlist(temp$Phat), ncol = n)
k <- temp$k

df_test<- data.frame(
  Graph <- I(list(G)),
  node = n,
  k = k,
  m = I(list(m1)),
  Phat = I(list(Phat)),
  fit = T
)

result <- generate_SBM(df_test)
test <- return_bootstrap_statistic(n, result, "SBM")
if ( (test < qtw(1 - (alpha / 2))) && (test > qtw(alpha / 2)) ) {reject[1] = 0}

# SBM 
m.file <- read.csv(paste("SBM_m_n", n, ".csv", sep = ""))
m <- unlist(m.file[i, ])
k = 3
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
  B = I(list(B)),
  model = "SBM"
)

temp <- SBM.mle(df_true)
G <- matrix(unlist(temp$G), ncol = n)
m1 <- unlist(temp$m)
Phat <- matrix(unlist(temp$Phat), ncol = n)

df_test<- data.frame(
  Graph <- I(list(G)),
  node = n,
  k = 3,
  m = I(list(m1)),
  Phat = I(list(Phat)),
  fit = T
)

result <- generate_SBM(df_test)
test <- return_bootstrap_statistic(n, result, "SBM")
if ( (test < qtw(1 - (alpha / 2))) && (test > qtw(alpha / 2)) ) {reject[2] = 0}

# Beta
file <- read.csv(paste("beta_n", n, ".csv", sep = ""))
vec <- t(t(file[i, ]))

df_true = data.frame(
  node = n,
  beta <- I(list(vec)),
  link <- "expit",
  model = "Beta"
)

temp <- SBM.mle(df_true)
G <- matrix(unlist(temp$G), ncol = n)
m1 <- unlist(temp$m)
Phat <- matrix(unlist(temp$Phat), ncol = n)
k <- temp$k
df_test<- data.frame(
  Graph <- I(list(G)),
  node = n,
  k = k,
  m = I(list(m1)),
  Phat = I(list(Phat)),
  fit = T
)

result <- generate_SBM(df_test)
test <- return_bootstrap_statistic(n, result, "SBM")
if ( (test < qtw(1 - (alpha / 2))) && (test > qtw(alpha / 2)) ) {reject[3] = 0}



# Plug in your csv file name here
