# Run simulations using given functions. Do a bootstrap correction if possible
library(RMTstat)
run_simulations <- function (node, df_true, df_test, true_model, test_model, alpha) {
  
  source(paste("Generate_", true_model, ".R", sep = "")) 
  source(paste("Generate_", test_model, ".R", sep = ""))
  G <- eval(parse(text = paste("generate_", true_model, "(df_true)", sep = "")))
  df_test <- cbind(df_test, Graph = I(list(G)))
  result <- eval(parse(text = paste("generate_", test_model, "(df_test)", sep = "")))
  source("Return_Bootstrap_Statistic.R")
  model = test_model
  test <- return_bootstrap_statistic(df_test$node, result, model)
  reject = 1
  if ( (test < qtw(1 - (alpha / 2))) && (test > qtw(alpha / 2)) ) {reject = 0}
  return (reject)
}

