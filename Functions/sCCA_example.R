#' An example running the functions presented in the repository
#' 
#' These functions were used in producing the results for the manuscript called
#' Integrating multi-OMICS data through sparse canonical correlation analysis for the prediction of complex traits: a comparison study, 
#' by Rodosthenous, Shahrezaei and Evangelou, in Bioinformatics
#' 

# Directory # 
#setwd("")
setwd("C:/Users/theod/OneDrive - Imperial College London/Documents/Imperial/PhD/R")

# Libraries # 
library(PMA)
source("ConvCCA.R")
source("RelPMDCCA.R")

# Simulate two datasets #
# Follow the example in PMC -- CCA function

# A simple simulated example
set.seed(1234)
u <- matrix(c(rep(1,25),rep(0,75)),ncol=1)
v1 <- matrix(c(rep(1,50),rep(0,450)),ncol=1)
v2 <- matrix(c(rep(0,50),rep(1,50),rep(0,900)),ncol=1)
x <- u%*%t(v1) + matrix(rnorm(100*500),ncol=500)
z <- u%*%t(v2) + matrix(rnorm(100*1000),ncol=1000)

# Run PMDCCA #
# Follow example in the documentation of function CCA in package PMA
pmdCCA_ex <- CCA(x,z,typex="standard",typez="standard",K=1)

# Run ConvCCA #
convCCA_ex <- convCCA(X_1 = x,
                   X_2 = z,
                   tauW_1 = 0.1,
                   tauW_2 = 0.1, 
                   nIter = 100)

# Run RelPMDCCA # 
relPMDCCA_ex <- relPMDCCA(X_1 = x,
                          X_2 = z,
                          lambda = 10,
                          tauW_1 = 0.8,
                          tauW_2 = 0.8,
                          nIter = 100)

# Plots
par(mfrow = c(1,3))
plot(pmdCCA_ex$u)
plot(convCCA_ex$W_1)
plot(relPMDCCA_ex$W_1)
