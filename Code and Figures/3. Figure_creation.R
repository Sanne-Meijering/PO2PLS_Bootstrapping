#### Setup
source("Functions.R")
load("/Datasets/Joint_values.Rdata")
library(xtable)
K <- 1000 # Number of permutations
sign <- c(0.001, 0.005, 0.01, 0.05) # Significance levels
results_template <- data.frame(Sensitivity = rep(0,4), Specificity = rep(0,4), FDR = rep(0,4),
                               row.names=c("0.001", "0.005", "0.01", "0.05")) # Template for results


#### Figure of PO2PLS results of the original dataset (data not included in GitHub repository)
#load("/Datasets/PO2PLS_chr17.RData")
#pdf("/Figures/PO2PLS_res.pdf", height=6, width=6)
#plot(Psteps2000$parameters$W[,1], ylab = "Joint loading")
#dev.off()

#### Peak in specific condition (test only)
# plot(Psteps2000$parameters$Wo[,2])
#### SVD of peak (test only)
# load("/Datasets/Xy_chr17.RData")
# SVD <- svd(X, 0, 10)
# plot(SVD$v[,1])
# plot(SVD$v[,2])

#### Histogram of non-parametric bootstrap in joint and low noise scenario
load("/Datasets/nonparaboot_joint.Rdata") # Load dataset
N_genes <- nrow(orig_res$parameters$W) # Number of genes
p_val <- nonpara_pval(npboot_res, N_genes, K) # Calculate p-values

# Create histogram
pdf("/Figures/hist_nonpara_joint_lnoise.pdf", height=6, width=6)
hist(p_val, xlim=c(0,1), xlab = "P-value")
dev.off()

#### Parametric bootstrapping table and ROC-plot (high noise)
load("/Datasets/paraboot_joint.Rdata") # Load dataset with joint components
N_genes <- nrow(orig_res$parameters$W) # Number of genes
X <- rep(0,N_genes) # Labels
X[high_joint] <- 1 # Set labels of genes with joint componenent to 1
p_val <- para_pval(orig_res, boot_res, N_genes, K) # Calculate p-values
results <- results_template

# Select significant values and calculate relevant statistics
for(i in 1:4){
  sign_p <- which(p_val < sign[i]) # Significant values
  TP <- sum(high_joint %in% sign_p) # Number of True Positives
  FN <- length(high_joint) - TP # Number of False Negatives
  FP <- length(sign_p) - TP # Number of False Positives
  TN <- N_genes - (TP+FN+FP) # Number of True Negatives
  
  # Calculate statistics
  results$Sensitivity[i] <- TP/(TP+FN)
  results$Specificity[i] <- TN/(TN+FP)
  results$FDR[i] <- FP/(TP+FP)
}

# Create ROC plot
pdf("/Figures/ROC_para_hnoise.pdf", height=6, width=6)
plot_roc(p_val, X)
dev.off()

load("/Datasets/paraboot_njoint.Rdata") # Load dataset without joint components
p_val <- para_pval(orig_res, boot_res, N_genes, K) # Calculate p-values

# Select significant values and calculate specificity
for(i in 1:4){
  sign_p <- which(p_val < sign[i])
  FP <- length(sign_p) # Number of False Positives
  TN <- N_genes - FP # Number of True Negatives
  results$"Specificity (no joint)"[i] <- TN/(TN+FP)
}

xtable(results, type="latex") # Create latex table

##### Permutation testing table and ROC-plot (high noise)
load("/Datasets/perm_joint.Rdata") # Load dataset with joint component
N_genes <- nrow(orig_res$parameters$W) # Number of genes
X <- rep(0,N_genes) # Labels
X[high_joint] <- 1 # Set labels of genes with joint componenent to 1
p_val <- ratio_calc(orig_res, perm_res, N_genes, K) # Calculate p-values
results <- results_template

# Select significant values and calculate relevant statistics
for(i in 1:4){
  sign_p <- which(p_val < sign[i]) # Significant values
  TP <- sum(high_joint %in% sign_p) # Number of True Positives
  FN <- length(high_joint) - TP # Number of False Negatives
  FP <- length(sign_p) - TP # Number of False Positives
  TN <- N_genes - (TP+FN+FP) # Number of True Negatives
  
  results$Sensitivity[i] <- TP/(TP+FN)
  results$Specificity[i] <- TN/(TN+FP)
  results$FDR[i] <- FP/(TP+FP)
}

# Create ROC plot
pdf("/Figures/ROC_perm_hnoise.pdf", height=6, width=6)
plot_roc(p_val, X)
dev.off()

load("/Datasets/perm_njoint.Rdata") # Load dataset without joint component
p_val <- ratio_calc(orig_res, perm_res, N_genes, K) # Calculate p-values

# Select significant values and calculate specificity
for(i in 1:4){
  sign_p <- which(p_val < sign[i]) # Significant values
  FP <- length(sign_p) # Number of False Positives
  TN <- N_genes - FP # Number of True Negatives
  results$"Specificity (no joint)"[i] <- TN/(TN+FP)
}

xtable(results, type="latex") # Create table

#### Parametric bootstrapping table and ROC-plot (low noise)
load("/Datasets/paraboot_joint_lnoise.Rdata") # Load dataset with joint component
N_genes <- nrow(orig_res$parameters$W) # Number of genes
X <- rep(0,N_genes) # Labels
X[high_joint] <- 1 # Set labels of genes with joint componenent to 1
p_val <- para_pval(orig_res, boot_res, N_genes, K) # Calculate p-values
results <- results_template

# Select significant values and calculate relevant statistics
for(i in 1:4){
  sign_p <- which(p_val < sign[i]) # Significant values
  TP <- sum(high_joint %in% sign_p) # Number of True Positives
  FN <- length(high_joint) - TP # Number of False Negatives
  FP <- length(sign_p) - TP # Number of False Positives
  TN <- N_genes - (TP+FN+FP) # Number of True Negatives
  
  results$Sensitivity[i] <- TP/(TP+FN)
  results$Specificity[i] <- TN/(TN+FP)
  results$FDR[i] <- FP/(TP+FP)
}

# Create ROC plot
pdf("/Figures/ROC_para_lnoise.pdf", height=6, width=6)
plot_roc(p_val, X)
dev.off()

load("/Datasets/paraboot_njoint_lnoise.Rdata") # Load dataset without joint component
p_val <- para_pval(orig_res, boot_res, N_genes, K) # Calculate p-values

# Select significant values and calculate specificity
for(i in 1:4){
  sign_p <- which(p_val < sign[i]) # Significant values
  FP <- length(sign_p) # Number of False Positives
  TN <- N_genes - FP # Number of True Negatives
  results$"Specificity (no joint)"[i] <- TN/(TN+FP)
}

xtable(results, type="latex") # Create table

##### Table of Permutation, joint & Low noise
load("/Datasets/perm_joint_lnoise.Rdata") # Load dataset
N_genes <- nrow(orig_res$parameters$W) # Number of genes
X <- rep(0,N_genes) # Labels
X[high_joint] <- 1 # Set labels of genes with joint componenent to 1
p_val <- ratio_calc(orig_res, perm_res, N_genes, K) # Calculate the p-values
results <- results_template

# Select significant values and calculate relevant statistics
for(i in 1:4){
  sign_p <- which(p_val < sign[i])
  TP <- sum(high_joint %in% sign_p) # Number of True Positives
  FN <- length(high_joint) - TP # Number of False Negatives
  FP <- length(sign_p) - TP # Number of False Positives
  TN <- N_genes - (TP+FN+FP) # Number of True Negatives
  
  results$Sensitivity[i] <- TP/(TP+FN)
  results$Specificity[i] <- TN/(TN+FP)
  results$FDR[i] <- FP/(TP+FP)
}

# Create ROC plot
pdf("/Figures/ROC_perm_lnoise.pdf", height=6, width=6)
plot_roc(p_val, X)
dev.off()

load("/Datasets/perm_njoint_lnoise.Rdata") # Load dataset without joint component
p_val <- ratio_calc(orig_res, perm_res, N_genes, K) # Calculate p-values

# Select significant values and calculate specificity
for(i in 1:4){
  sign_p <- which(p_val < sign[i])
  FP <- length(sign_p) # Number of False Positives
  TN <- N_genes - FP # Number of True Negatives
  results$"Specificity (no joint)"[i] <- TN/(TN+FP)
}

xtable(results, type="latex") # Create latex table