#### Setup
source("Functions.R")
load("Datasets.Rdata")
steps <- 100 # Number of steps for PO2PLS
K <- 1000 # Number of bootstrapping iterations
r <- 1 # Number of joint components
rx <- 2 # Number of specific components of X
ry <- 4 # Number of specific components of Y
n.cores <- 1 # Number of cores for parallel bootstrapping

#### Scenario with joint component, high noise
# Set dataset
set.seed(26452)
data <- data_joint
N <- nrow(data$X) # Sample size

# Run simulation
orig_res <- PO2PLS(data$X,data$Y,r,rx,ry, steps=steps) # Run PO2PLS on original data
npboot_res <- nonparaboot(orig_res, data$X, data$Y, K, steps, boot.cores=n.cores) # Do non-parametric bootstrapping
save(orig_res, npboot_res, file="/Datasets/nonparaboot_joint.Rdata") # Store results

#### Scenario without joint component, high noise
# Set dataset
set.seed(26452)
data <- data_njoint
N <- nrow(data$X) # Sample size

# Run simulation
orig_res <- PO2PLS(data$X,data$Y,r,rx,ry, steps=steps) # Run PO2PLS on original data
npboot_res <- nonparaboot(orig_res, data$X, data$Y, K, steps, boot.cores=n.cores) # Do non-parametric bootstrapping
save(orig_res, npboot_res, file="/Datasets/nonparaboot_njoint.Rdata") # Store results

#### Scenario with joint component, low noise
# Set dataset
set.seed(26452)
data <- data_joint_lnoise
N <- nrow(data$X) # Sample size

# Run simulation
orig_res <- PO2PLS(data$X,data$Y,r,rx,ry, steps=steps) # Run PO2PLS on original data
npboot_res <- nonparaboot(orig_res, data$X, data$Y, K, steps, boot.cores=n.cores) # Do non-parametric bootstrapping
save(orig_res, npboot_res, file="/Datasets/nonparaboot_joint_lnoise.Rdata") # Store results

#### Scenario without joint component, low noise
# Set dataset
set.seed(26452)
data <- data_njoint_lnoise
N <- nrow(data$X) # Sample size

# Run simulation
orig_res <- PO2PLS(data$X,data$Y,r,rx,ry, steps=steps) # Run PO2PLS on original data
npboot_res <- nonparaboot(orig_res, data$X, data$Y, K, steps, boot.cores=n.cores) # Do non-parametric bootstrapping
save(orig_res, npboot_res, file="/Datasets/nonparaboot_njoint_lnoise.Rdata") # Store results