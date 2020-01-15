# Generate datasets

# Setup
library(OmicsPLS)
library(PO2PLS)
library(parallel)
library(tidyverse)
library(nipals)
set.seed(42342)

# Load dataset (smallest for now)
setwd("C:/Users/rikul/Documents/Thesis_data")
load("Xy_chr21.Rdata")
res <- nipals(X, 20, verbose=T)
load("components20_chr21.Rdata")
plot(res$eig, type="b", xlab="Component", ylab="Eigenvalue", ylim=c(0,230)) 
abline(h=1, col="red")
# The elbow seems to be at co

# Set variables
N <- 100 # Sample size
p <- 100 # Number of X variables
q <- 7 # Number of Y variables
K <- 50 # Number of bootstrapping iterations
r <- 1 # Number of joint components
rx <- 3 # Number of specific components of X
ry <- 4 # Number of specific components of Y

# Generate parameters given the set variables
params <- generate_params(p, q, r, rx, ry, type = "unit")

# Change parameters for different datasets
# Dataset 1: No joint components, but 