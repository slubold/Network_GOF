
library(ergm)
library(sna)
library(RMTstat)
source("Generate_ERGM.R")
df = data.frame(
  node = 500,
  formula = "edges + triangle",
  param = I(list(c(-1/2, -1/3)))
)
G <- generate_ERGM(df)

df1 <- data.frame(
  node = 500,
  formula = "edges",
  Graph = I(list(G)),
  fit = T
)
result <- generate_ERGM(df1)
result$MLE

# Save to a csv file
n <- 500
source("Return_Bootstrap_Statistic.R")
test <- return_bootstrap_statistic(n, result, model = "ERGM")
test
