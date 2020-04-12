###############################################################################
########### Simulations ###########
###############################################################################

library(MASS)
library(data.table)
library(corpcor)

make.positive.definite <- function (m, tol)
{ # From caret library
    if (!is.matrix(m))
        m = as.matrix(m)
    d = dim(m)[1]
    if (dim(m)[2] != d)
        stop("Input matrix is not square!")
    es = eigen(m, symmetric = TRUE)
    esv = es$values
    if (missing(tol))
        tol = d * max(abs(esv)) * .Machine$double.eps
    delta = 2 * tol
    tau = pmax(0, delta - esv)
    dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
    return(m + dm)
}

# Name covariance measures
## 1. identity
covarianceIdentity <- function(p){
  Sigma <- diag(p)
  return(Sigma)
}
## 2. toeplitz
covarianceToeplitz <- function(p){
  Sigma <- matrix(0, nrow = p, ncol =p)
  for (i in 1:p){
    for (j in 1:p){
      power <- abs(i-j)
      Sigma[i,j] = 0.9^power
    }
  }
  return(Sigma)
}
## 3. sparse-inverse
covarianceSparse <- function(p){
  # Define Omega
  Sigma <- matrix(0, nrow = p, ncol = p)
  Omega <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p){
    for (j in 1:p){
      power <- abs(i-j)
      if (power == 0){
        Omega[i,j] = 1
      } else if (power == 1){
        Omega[i,j] = 0.5
      } else if (power == 2){
        Omega[i,j] = 0.4
      } else {
        Omega[i,j] = 0
      }
    }
  }
  Sigma0 <- solve(Omega)
  for (i in 1:p){
    for (j in 1:p){
      denominator <- sqrt(Sigma0[i,i]*Sigma0[j,j])
      Sigma[i,j] = Sigma0[i,j]/denominator
    }
  }
  return(Sigma)
}

simulationCovariance <- function(type = "identity", n, p, q, rho = 0.9, uTrue, vTrue){
  # type can take the following: "identity", "toeplitz", "sparse", "spiked"
  # Generate the covariance matrix:
  if (type == "identity"){
    SigmaX <- covarianceIdentity(p)
    SigmaY <- covarianceIdentity(q)
  } else if (type == "toeplitz"){
    SigmaX <- covarianceToeplitz(p)
    SigmaY <- covarianceToeplitz(q)
  } else if (type == "sparse"){
    SigmaX <- covarianceSparse(p)
    SigmaY <- covarianceSparse(q)
  }
  SigmaXY <- rho*(SigmaX%*%uTrue%*%t(vTrue)%*%SigmaY)
  SigmaYX <- t(SigmaXY)
  ## Generate the data-sets
  Sigma1 <- cbind(SigmaX, SigmaXY)
  Sigma2 <- cbind(SigmaYX, SigmaY)
  Sigma <- rbind(Sigma1, Sigma2)
  Sigma <- make.positive.definite(Sigma)
  XY <- mvrnorm(n = n, mu = rep(0,p+q), Sigma = Sigma)
  X <- XY[,1:p]
  Y <- XY[,-(1:p)]
  return(list("X" = X, "Y" = Y))
}
## Obtain the initial guess:
initialGuess <- function(X, Y, thres){
  # t: thresholding
  # Define matrix:
  tempMatrix <- t(X)%*%Y
  S <- matrix(0,nrow = nrow(tempMatrix), ncol = ncol(tempMatrix))
  ## Soft-thresholding, element-wise
  for (i in 1:nrow(tempMatrix)){
    for (j in 1:ncol(tempMatrix)){
      temp <- tempMatrix[i,j]
      if (abs(temp) > thres){
        S[i,j] = sign(temp)*(abs(temp) - thres)
      } else {
        S[i,j] = 0
      }
    }
  }
  svdS <- svd(S)
  Uhat <- svdS$u
  Vhat <- svdS$v
  Utilde <- normalizeForInitial(Uhat, X)
  Vtilde <- normalizeForInitial(Vhat, Y)
  # Calculate Dtilde
  Dtilde = t(Utilde)%*%t(X)%*%Y%*%Vtilde
  k = which.max(diag(Dtilde))
  uInitial <- Utilde[,k]
  vInitial <- Vtilde[,k]
  return(list("uInit" = uInitial, "vInit" = vInitial))
}
## Normalization as suggested by Suo et al. (RelPMDCCA)
normalizeForInitial <- function(u, X){
  newU <- matrix(0, nrow = nrow(u), ncol = ncol(u))
  newColU <- numeric(length = nrow(u))
  for (i in 1:ncol(u)){
    colU <- as.vector(u[,i])
    denominator <- sqrt(t(colU)%*%(t(X)%*%X)%*%colU)
    for (j in 1:length(colU)){
      newColU[j] = colU[j]/denominator
    }
    newU[,i] <- newColU
  }
  return(newU)
}

simulationSingleLatent <-  function(n, p, p.cc, q, q.cc, sd.mu, sd.epsilon){
  # generate alpha - coefficients of mu for X
  alpha <- runif(p.cc)
  # generate beta - coefficients of mu for Y
  beta <- runif(q.cc)
  # Normalize
  alpha <- alpha/sum(alpha)
  beta <- beta/sum(beta)
  # generate mu
  mu <- rnorm(n = n, mean = 0, sd = sd.mu) # Simulate mu
  # Initiate
  X <- matrix(nrow = n, ncol = p)
  Y <- matrix(nrow = n, ncol = q)
  ## Compute X:
  for (i in 1:nrow(X)){
    for (j in 1:ncol(X)){
      if (j <= p.cc){
        X[i,j] <- alpha[j]*mu[i] + rnorm(n=1, mean = 0, sd = sd.epsilon)
      } else {
        X[i,j] <- rnorm(n=1, mean = 0, sd = sd.epsilon)
      }
    }
  }
  ## Compute Y:
  for (i in 1:nrow(Y)){
    for (j in 1:ncol(Y)){
      if (j <= q.cc){
        Y[i,j] <- beta[j]*mu[i] + rnorm(n=1, mean = 0, sd = sd.epsilon)
      } else {
        Y[i,j] <- rnorm(n=1, mean = 0, sd = sd.epsilon)
      }
    }
  }
  # Define mu.x and mu.y
  muX <- rowSums(X[,1:p.cc])
  muY <- rowSums(Y[,1:q.cc])
  ## Move on to computing the Covariance matrix in sets X and Y
  ## Note that it has a block structure:
  SigmaXY <- matrix(nrow = p, ncol = q)
  for (i in 1:nrow(SigmaXY)){
    for (j in 1:ncol(SigmaXY)){
      if ( (i <= p.cc) & (j <= q.cc) ){
        SigmaXY[i,j] <- alpha[i]*beta[j]*(sd.mu^2)
      } else {
        SigmaXY[i,j] <- 0
      }
    }
  }
  ## If true u, v and rho are computed through SVD:
  svdSigmaXY <- svd(SigmaXY)
  uTrue <- svdSigmaXY$u[,1]
  vTrue <- svdSigmaXY$v[,1]
  rhoTrue <- svdSigmaXY$d[1]
  return(list("X" = X, "Y" = Y,
              "uTrue" = uTrue, "vTrue" = vTrue, "rhoTrue" = rhoTrue))
}

simulationSingleLatentMultiple <-  function(n, pList, p.ccList, sd.mu, sd.epsilon){
  # number of data-sets
  length.list <- length(pList)
  if (length(pList) != length(p.ccList)){
    break("Please provide information on all data-sets, with the length of features being the same as the length of cross-correlated features")
  }
  alpha <- vector("list", length = length.list)
  for (i in 1:length.list){
    alpha[[i]] <- runif(p.ccList[[i]])
    alpha[[i]] <- alpha[[i]]/sum(alpha[[i]]) # sum to 1
  }
  # Generate mu
  rnorm(n = n, mean = 0, sd = sd.mu)
  # Initiate
  X <- vector("list", length = length.list)
  for (i in 1:length.list){
    X[[i]] <- matrix(nrow = n, ncol = p[[i]])
  }
  ## Compute Xlist:
  for (l in 1:length.list){
    for (i in 1:nrow(X[[l]])){
      for (j in 1:ncol(X[[l]])){
        if (j <= p.ccList[[l]]){
          X[[l]][i,j] <- alpha[[l]][j]*mu[i] + rnorm(n = 1, mean = 0, sd = sd.epsilon)
        } else {
          X[[l]][i,j] <- rnorm(n = 1, mean = 0, sd = sd.epsilon)
        }
      }
    }
  }
  # Define mu.xList
  mu.xList <- vector("list", length = length.list)
  for (i in 1:length.list){
    mu.xList[[i]] <- rowSums(X[[i]][,1:p.ccList[[i]]])
  }
  return(list("X" = X))
}

runSimulation <- function(method = "singleLatent", parameters){
  if (method == "singleLatent"){
    ## Take parameters:
    n = parameters[1]
    p = parameters[2]
    p.cc = parameters[3]
    q = parameters[4]
    q.cc = parameters[5]
    sd.mu = parameters[6]
    sd.epsilon = parameters[7]
    # Simulation
    simulPark <- simulationParkhomenko(n, p, p.cc, q, q.cc, sd.mu, sd.epsilon)
    X <- simulPark$X
    Y <- simulPark$Y
    uTrue <- simulPark$uTrue
    vTrue <- simulPark$vTrue
    rhoTrue <- simulPark$rhoTrue
  } else if (method == "covarianceBased"){
    type = parameters[1]
    n = as.numeric(parameters[2])
    p = as.numeric(parameters[3])
    q = as.numeric(parameters[4])
    rhoTrue = as.numeric(parameters[5])
    numNZu = as.numeric(parameters[6])
    numNZv = as.numeric(parameters[7])
    uTrue <- c(runif(numNZu), rep(0, (p-numNZu)))
    vTrue <- c(runif(numNZv), rep(0, (q-numNZv)))
    simulSuo <- simulationSuo(type, n, p, q, rhoTrue, uTrue, vTrue)
    X <- simulSuo$X
    Y <- simulSuo$Y
  }  else {
    stop("Please provide a valid simulation. Choose from singleLatent, covarianceBased")
  }
  return(list("X" = X, "Y" = Y, "uTrue" = uTrue,
              "vTrue" = vTrue, "rhoTrue" = rhoTrue))
}
