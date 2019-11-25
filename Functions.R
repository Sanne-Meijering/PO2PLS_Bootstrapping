# Cross-validation for PO2PLS
max_logl_PO2PLS <- function(X, y, r, rx, ry, steps=100, init_param="unit"){
  dim_names <- list(r=1:r, rx=0:rx, ry=0:ry)
  total <- r*(rx+1)*(ry+1)
  likelihood <- array(rep(0, total), dim=c(r, rx+1, ry+1), dimnames = dim_names)
  for(i in 1:r){
    for(j in 0:rx){
      for(k in 0:ry){
        if((r + rx <= ncol(X)) & (r + ry <= ncol(y))){
          res <- PO2PLS(X, y, r, rx, ry, steps = steps, init_param = "unit", verbose=F)
          likelihood[i, j+1, k+1] <- res$meta_data$loglikelihood$last_val
        }
      }
    }
  }
  maximum = which(likelihood==min(likelihood), arr.ind=T)
  maximum[-1] <- maximum[-1] - 1
  
  return(list(maximum = maximum, likelihood_array = likelihood))
}