# X = X-dataset
# y = y-dataset
# r, rx, ry = respectively number of joint components, x-specific components & y-specific components
# K = number of bootstrapping iterations
# type = type of bootstrapping
nonparaboot <- function(X, Y, r, rx, ry, K){
  # Create list for results
  bootstrap_res <- list()
  
  # Set sizes of samples and dimensions of dataset
  N <- nrow(X)
  NX <- ncol(X)
  NY <- ncol(Y)
  
  for(i in 1:K){
    #
    bootsamp <- sample(1:N, 1:N, replace=T)
    sampX <- X[bootsamp,]
    sampY <- Y[bootsamp,]
    res <- PO2PLS(sampX, sampY, r, rx, ry, steps=2000)
  }
  
}

# PO2PLS_param = parameters from PO2PLS results
# K = number of bootstrap iterations
# N = sample size of bootstraps
paraboot <- function(PO2PLS_param, K, N){
  rnorm(2000)
}