#### Setup
source("Functions.R")
load("/Datasets/Datasets.Rdata")
# Set variables
K <- 1000 # Number of permutation iterations
r <- 1 # Number of joint components
rx <- 2 # Number of specific components of X
ry <- 4 # Number of specific components of Y
n.cores <- 1 # Number of cores for parallel permutation testing

#### Scenario with joint component, high noise
set.seed(26452)
data <- data_joint
X <- data$X
y <- data$Y

# Run simulation
orig_res <- PO2PLS(X,y,r,rx,ry, steps=100) # Run PO2PLS on original data
perm_res <- PO2PLSpermute_multi(X,y,r,rx,ry,K,perm.cores=n.cores) # Do permutation testing
save(orig_res, perm_res, file="/Datasets/perm_joint.Rdata") # Store results

#### Scenario without joint component, high noise
set.seed(26452)
data <- data_njoint
X <- data$X
y <- data$Y

# Run simulation
orig_res <- PO2PLS(X,y,r,rx,ry, steps=100) # Run PO2PLS on original data
perm_res <- PO2PLSpermute_multi(X,y,r,rx,ry,K,perm.cores=n.cores) # Do permutation testing
save(orig_res, perm_res, file="/Datasets/perm_njoint.Rdata") # Store results

#### Scenario with joint component, low noise
set.seed(26452)
data <- data_joint_lnoise
X <- data$X
y <- data$Y

# Run simulation
orig_res <- PO2PLS(X,y,r,rx,ry, steps=100) # Run PO2PLS on original data
perm_res <- PO2PLSpermute_multi(X,y,r,rx,ry,K,perm.cores=n.cores) # Do permutation testing
save(orig_res, perm_res, file="/Datasets/perm_joint_lnoise.Rdata") # Store results

#### Scenario without joint component, low noise
set.seed(26452)
data <- data_njoint_lnoise
X <- data$X
y <- data$Y

# Run simulation
orig_res <- PO2PLS(X,y,r,rx,ry, steps=100) # Run PO2PLS on original data
perm_res <- PO2PLSpermute_multi(X,y,r,rx,ry,K,perm.cores=n.cores) # Do permutation testing
save(orig_res, perm_res, file="/Datasets/perm_njoint_lnoise.Rdata") # Store results