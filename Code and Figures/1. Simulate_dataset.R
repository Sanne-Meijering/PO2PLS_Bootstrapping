#### Setup
library(OmicsPLS)
library(PO2PLS)
library(dplyr)
set.seed(31*12*2016)
# Load the chromosome 17 data
load("/Datasets/Xy_chr17.RData")

#### Obtain parameters
# Select the optimal number of components
# Crossvalidation with max. 3 components
crossval_o2m_adjR2(X, y, 1:3, 0:3, 0:3, nr_folds=5)
# n=1, nx=2 and ny=3 for top result, ny=3 for all 3 top results
crossval_o2m_adjR2(X, y, 1:2, 2:4, 2:4, nr_folds=5)
# nx=2:4 and ny=4, repeat with higher max ny
crossval_o2m_adjR2(X, y, 1, 2:5, 3:5, nr_folds=5)
# n=1 nx=2 & ny=4 for best results, same results as for previous test

#### Remove peak area (testing only)
# X <- X[,c(1:3100, 3700:6702)]
# save(Psteps2000, "Psteps2000", file="/Datasets/PO2PLS_ch17_new.Rdata")
# plot(Psteps2000$parameters$W[,1])

# Run PO2PLS
Psteps2000 <- PO2PLS(X, y, 1, 2, 4, steps=2000)
save(Psteps2000, "Psteps2000", file="/Datasets/PO2PLS_ch17.Rdata")

#### Simulation of datasets
# Load data and set seed
load("/Datasets/PO2PLS_chr17.RData")
set.seed(3127391)

# Set parameters
N <- 1000 # Number of participants
params <- Psteps2000$parameters # Store parameters
r <- ncol(params$W) # Number of joint components
rx <- ncol(params$Wo) # Number of specific components of X
ry <- ncol(params$Co) # Number of specific components of Y
N_genes <- nrow(params$W) # Number of genes
N_hist <- nrow(params$C) # Number of histological traits

# There is a major peak in the joint variables between 3100 and 3700. 
# These joint variables were remove to avoid them from affecting results
params$W <- as.matrix(params$W[c(0:3100, 3700:N_genes),], ncol=r)
params$Wo <- as.matrix(params$Wo[c(0:3100, 3700:N_genes),], ncol=r)
# Reset the number of genes
N_genes <- nrow(params$W)

# Remove joint components (condition without joint component)
params_njoint <- params
params_njoint$W <- as.matrix(rep(0, N_genes*r), ncol=r)
params_njoint$C <- as.matrix(rep(0, N_hist*r), ncol=r)

# Remove all but the highest 200 joint components (condition with joint component)
high_joint <- sort(params$W, decreasing = T, index.return=T)$ix[1:200]
save(high_joint, file="Joint_values.Rdata")
no_joint <- 1:N_genes
no_joint <- no_joint[-high_joint]
params_joint <- params
params_joint$W[no_joint,] <- as.matrix(rep(0, (N_genes-200)*r), ncol=r)

# Split up in two different noise levels
# A signal-to-noise ratio of 4:1 is wanted, this means sigE should be:
sig2E <- (PO2PLS::tr(params$SigT) + PO2PLS::tr(params$SigTo))/4/N_genes

# Create lnoise datasets
params_njoint_lnoise <- params_njoint
params_njoint_lnoise$sig2E <- sig2E
params_joint_lnoise <- params_joint
params_joint_lnoise$sig2E <- sig2E

# Generate datasets
data_njoint <- PO2PLS:::generate_data(N, params_njoint)
data_joint <- PO2PLS:::generate_data(N, params_joint)
data_njoint_lnoise <- PO2PLS:::generate_data(N, params_njoint_lnoise)
data_joint_lnoise <- PO2PLS:::generate_data(N, params_joint_lnoise)

# Save datasets
save(data_joint, data_njoint, data_joint_lnoise, data_njoint_lnoise, 
     list=c("data_joint", "data_njoint", "data_joint_lnoise", 
            "data_njoint_lnoise"), file="/Datasets/Datasets.Rdata")