# Find number of components with highest Maximum Likelihood for PO2PLS 
max_logl_PO2PLS <- function(X, y, r, rx, ry, steps=100, init_param="unit"){
  X <- as.matrix(X)
  y <- as.matrix(y)
  ny <- ncol(y)
  nX <- ncol(X)
  dim_names <- list(r=1:r, rx=0:rx, ry=0:ry)
  total <- r*(rx+1)*(ry+1)
  likelihood <- array(rep(0, total), dim=c(r, rx+1, ry+1), dimnames = dim_names)
  for(i in 1:r){
    for(j in 0:rx){
      for(k in 0:ry){
        if((i + j <= min(nX, ny)) & (i + k <= min(nX, ny))){
          res <- PO2PLS(X, y, i, j, k, steps = steps, init_param = "o2m")
          likelihood[i, j+1, k+1] <- res$logl[[1]]
        }
      }
    }
  }
  maximum = which(likelihood==min(likelihood), arr.ind=T)
  maximum[-1] <- maximum[-1] - 1
  
  return(list(maximum = maximum, likelihood_array = likelihood))
}

# Cross-validation for PO2PLS
cross_val_PO2PLS <- function(X, y, r, rx, ry, steps=100, init_param="o2m", n_fold=1){
  # Divide the dataset into n_fold groups
  fold_size <- floor(nrow(X)/n_fold) # Determine fold size
  folds <- list() # Create list for samples
  X_remaining <- 1:nrow(X) # Make list samples without fold
  
  for(i in n_fold - 1){
    # Sample fold
    folds[[i]] <- sample(X_remaining, fold_size, replace=F)
    # Remove sampled values 
    X_remaining <- X_remaining[-match(folds[[i]], X_remaining)]
  }
  # Last fold is what is remaining in X_remaining
  folds[[n_fold]] <- X_remaining
  
  # Create an array for the likelihoods
  dim_names <- list(r=1:r, rx=0:rx, ry=0:ry)
  total <- r*(rx+1)*(ry+1)
  likelihood <- array(rep(0, total), dim=c(r, rx+1, ry+1), dimnames = dim_names)
  
  # Go through all values of r, rx and ry
  for(i in 1:r){
    for(j in 0:rx){
      for(k in 0:ry){
        # Skip all values where the number of joint+specific components exceeds the number of variables.
        if((r + rx <= ncol(X)) & (r + ry <= ncol(y))){
          for(i in 1:n_fold){
            train_set <- folds[[]]
            res <- PO2PLS(folds[[]], y, r, rx, ry, steps = steps, init_param = "o2m", verbose=F)
            
          }
          
          likelihood[i, j+1, k+1] <- res$meta_data$loglikelihood$last_val
        }
      }
    }
  }
  maximum = which(likelihood==min(likelihood), arr.ind=T)
  maximum[-1] <- maximum[-1] - 1
  
  return(list(maximum = maximum, likelihood_array = likelihood))
}