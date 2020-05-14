library(parallel)
library(PO2PLS)
library(OmicsPLS)

#### A function for non-parametric bootstrapping with PO2PLS
# fit = PO2PLS fit
# X = X-dataset, y = y-dataset
# K = number of bootstrapping iterations
# steps = number of PO2PLS steps per bootstrapping iteration
# boot.cores = number of cores used
nonparaboot <- function(fit, X, Y, K, steps=100, boot.cores=1){
  # Set number of samples
  N <- nrow(X)
  # Set number of components
  r <- ncol(fit$parameters$W)
  rx <- ncol(fit$parameters$Wo)*sign(ssq(fit$parameters$Wo))
  ry <- ncol(fit$parameters$Co)*sign(ssq(fit$parameters$Co))
  # Run bootstrap on multiple cores
  boot.par <- mclapply(1:K,FUN = function(niks){
    # Sample from dataset with replacement
    perm.indx <- sample(N, replace = TRUE)
    # Fit with scaled versions of the permutation
    PO2PLS(scale(X[perm.indx, ]), scale(Y[perm.indx, ]), r, rx, ry, steps=steps)
  }, mc.cores = boot.cores)
}

#### A function for parametric bootstrapping with PO2PLS
# fit = PO2PLS fit
# K = number of bootstrap iterations
# N = sample size of bootstraps
# steps = number of PO2PLS steps per bootstrapping iteration
# boot.cores = number of cores used
paraboot <- function(fit, K, N, steps=100, boot.cores=1){
  # Set number of components
  r <- ncol(fit$parameters$W)
  rx <- ncol(fit$parameters$Wo)*sign(ssq(fit$parameters$Wo))
  ry <- ncol(fit$parameters$Co)*sign(ssq(fit$parameters$Co))
  N_genes <- nrow(fit$parameters$W) # Number of genetic PCs
  N_hist <- nrow(fit$parameters$C) # Number of histological traits
  # Set joint component to zero to create a null distribution
  fit$parameters$W <- matrix(rep(0, N_genes*r), ncol=r)
  fit$parameters$C <- matrix(rep(0, N_hist*r), ncol=r)
  # Run bootstrap on multiple cores
  boot.par <- mclapply(1:K,FUN = function(niks){
    # Create new dataset
    new_data <- PO2PLS:::generate_data(N, fit$parameters)
    # Fit with scaled versions of the bootstrapped dataset
    PO2PLS(scale(new_data$X), scale(new_data$Y), r, rx, ry, steps=steps)
  }, mc.cores = boot.cores)
}

#### Function for permutation testing with PO2PLS
# X = X dataset, Y = Y dataset
# r, rx, ry = Number of joint components and number of specific components for X and Y respectively
# K = Number of bootstrapping iterations
# perm.cores = Number of cores used
# steps = Number of PO2PLS steps per permutation iteration
PO2PLSpermute_multi <- function(X, Y, r, rx, ry, K, perm.cores=1, steps=100){  
  N <- nrow(X) # Number of genes
  # Run bootstrap on multiple cores
  boot.par <- mclapply(1:K,FUN = function(niks){
    # Create permutation
    perm.indx <- sample(N, replace = FALSE)
    # Fit with scaled versions of the permutation
    PO2PLS(scale(X), scale(y[perm.indx, ]), r, rx, ry, steps=steps)
  }, mc.cores = perm.cores)
}

#### Function for p-value calculation for permutation testing (absolute values used)
# orig_res = PO2PLS fit of original dataset
# perm_res = List of PO2PLS fits on permuted datasets
# N_genes = Number of genetic PCs
# K = Number of bootstrapping iterations
ratio_calc <- function(orig_res, perm_res, N_genes, K){
  pval <- rep(0,N_genes) # Create empty vector
  for(i in 1:N_genes){
    # Retrieve original value for measure of interest
    orig_val <- abs(orig_res$parameters$W[i,1]*orig_res$parameters$B[1,1])
    # Gather all permutation values and sort them
    perm_val <- rep(0,K)
    for(j in 1:K){
      perm_val[j]<- abs(perm_res[[j]]$parameters$W[i,1]*perm_res[[j]]$parameters$B[1,1])
    }
    # Calculate p-values:
    # Take the ratio of permuted values more extreme than the original value
    # Add half of the ratio of permuted values equal to the original value
    pval[i] <- mean(orig_val < perm_val) +.5*mean(perm_val == orig_val)
  }
  return(pval)
}

#### Function for p-value calculation for non-parametric bootstrapping (two-sided)
# boot_res = List of PO2PLS fits on bootstrapped datasets
# N_genes = Number of genetic PCs
# K = Number of bootstrapping iterations
# null.values = Expected value of measure of interest under the null hypothesis
nonpara_pval <- function(boot_res, N_genes, K, null.val=0){
  pval <- rep(0,N_genes) # Create empty vector
  for(i in 1:N_genes){
    # Gather all bootstrapping values and sort them
    boot_val <- rep(0,K)
    for(j in 1:K){
      boot_val[j]<- boot_res[[j]]$parameters$W[i,1]*boot_res[[j]]$parameters$B[1,1]
    }
    # Calculate p-values:
    # Take the ratio of bootstrapped values larger than the null value
    # Add half of the ratio of bootstrapped values equal to the null value
    half.pval <- mean(boot_val > null.val) +.5*mean(boot_val == null.val)
    # Convert to p-value by taking the lowest ratio (ratio or 1-ratio) and multiplying it by 2 
    pval[i] <- 2*min(c(half.pval,1-half.pval))
  }
  return(pval)
}

#### Function for p-value calculation for non-parametric bootstrapping (absolute values)
# boot_res = List of PO2PLS fits on bootstrapped datasets
# N_genes = Number of genetic PCs
# K = Number of bootstrapping iterations
# null.values = Expected value of measure of interest under the null hypothesis
nonpara_pval_abs <- function(boot_res, N_genes, K, null.val=0){
  pval <- rep(0,N_genes) # Create empty vector
  for(i in 1:N_genes){
    # Gather all bootstrapping values and sort them
    boot_val <- rep(0,K)
    for(j in 1:K){
      boot_val[j]<- abs(boot_res[[j]]$parameters$W[i,1]*boot_res[[j]]$parameters$B[1,1])
    }
    
    # Calculate p-values:
    # Take the ratio of bootstrapped values less extreme than the original value
    # Add half of the ratio of bootstrapped values equal to the original value
    pval[i] <- mean(boot_val < null.val) +.5*mean(boot_val == null.val)
  }
  return(pval)
}

#### Function for p-value calculation for parametric bootstrapping (absolute values)
# orig_res = PO2PLS fit of original dataset
# boot_res = List of PO2PLS fits on bootstrapped datasets
# N_genes = Number of genetic PCs
# K = Number of bootstrapping iterations
para_pval <- function(orig_res, boot_res, N_genes, K){
  pval <- rep(0,N_genes) # Create empty vector
  for(i in 1:N_genes){
    # Retrieve original value for measure of interest
    orig_val <- abs(orig_res$parameters$W[i,1]*orig_res$parameters$B[1,1])
    # Gather all bootstrapping values and sort them
    boot_val <- rep(0,K)
    for(j in 1:K){
      boot_val[j]<- abs(boot_res[[j]]$parameters$W[i,1]*boot_res[[j]]$parameters$B[1,1])
    }
    
    # Calculate p-values:
    # Take the ratio of bootstrapped values less extreme than the original value
    # Add half of the ratio of bootstrapped values equal to the original value
    pval[i] <- mean(orig_val < boot_val) +.5*mean(boot_val == orig_val)
  }
  return(pval)
}

#### Function for creating a ROC plot
# scores = p-values
# labels = labels: 0 = no joint component, 1 = joint component
plot_roc <- function(scores, labels){
  # Sort labels based on the p-values
  labels <- labels[order(scores)]
  # Create dataframe with the true positive rate (sensitivity) and false positive rate (1-specificity)
  roc <- data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels))
  # Create plot
  plot(roc$FPR, roc$TPR, type='l', xlab="1 - Specificity", ylab="Sensitivity")
}