###############################################################################
########### RelPMDCCA ###########
###############################################################################

# To raise a matrix x to a given n power
"%^%" <- function(x, n){
    with(eigen(x), vectors %*% (values^n * t(vectors)))
}

## Functions to update canonical vectors within convCCA algorithm
# penalty can take "SCAD", or "LASSO", or "ELASTIC-NET" all capital
updateW.relPMDCCA <- function(omega, mu, c, tau, tau_EN = NULL a, old, penalty = "LASSO", element_wise = TRUE){
    # tau_EN is the tuning parameter for the second element. It's used in the Elastic-Net
    # Step 1 -- Initialisation
    eta_1 <- -1/(2*(a-1))
    eta_2 <- (2*a*tau)/(2*(a-1))
    check <- omega + mu*c # Function criterion

    # Step 2 -- Implement algorithm via given penalty function
    names(check) <- 1:length(check)
    if ((penalty == "SCAD") %% (element_wise == FALSE)){
        # Split into five parts, as per the proximal criteria
        # First part
        condition_1 <- check[apply(check, 1, function(check){return((0 < check - mu*tau) && (check - mu*tau <= tau))})  ]
        updateW_1 <- condition_1 - mu*tau
        names(updateW_1) <- names(condition_1)
        # Second part
        condition_2 <- check[apply(check, 1, function(check){return((-tau <= check + mu*tau) && (check +mu*tau < 0))})]
        updateW_2 <- condition_2 + mu*tau
        names(updateW_2) <- names(condition_2)
        # Third part
        condition_3 <- check[apply(check, 1, function(check){return((tau  < (check - mu*eta_2)/(1+2*mu*eta_1)) && ((check - mu*eta_2)/(1+2*mu*eta_1) <= a*tau))})]
        updateW_3 <- (condition_3 - mu*eta_2)/(1+2*mu*eta_1)
        names(updateW_3) <- names(condition_3)
        # Fourth part
        condition_4 <- check[apply(check, 1, function(check){return((-a*tau  <= (check + mu*eta_2)/(1+2*mu*eta_1)) && ((check + mu*eta_2)/(1+2*mu*eta_1) <= -tau))})]
        updateW_4 <- (condition_4 + mu*eta_2)/(1+2*mu*eta_1)
        names(updateW_4) <- names(condition_4)
        # Fifth part
        condition_5 <- check[apply(check, 1, function(check){return(abs(check) > a*tau)})]
        updateW_5 <- condition_5
        names(updateW_5) <- names(condition_5)
        finalUpdate <- c(updateW_1, updateW_2, updateW_3, updateW_4, updateW_5)
        finalUpdate <- finalUpdate[names(check)] # Puts the data back to the original order
        finalUpdate <- as.numeric(finalUpdate) # Remove the names
    } else if (penalty = "SCAD") && (element_wise == TRUE)){
        # Split into five parts, as per the proximal criteria
        # First part
        condition_1 = (0 < check - mu*tau) && (check - mu*tau <= tau)
        condition_2 = (-tau <= check + mu*tau) && (check +mu*tau < 0)
        condition_3 = (tau  < (check - mu*eta_2)/(1+2*mu*eta_1)) && ((check - mu*eta_2)/(1+2*mu*eta_1) <= a*tau)
        condition_4 = (-a*tau  <= (check + mu*eta_2)/(1+2*mu*eta_1)) && ((check + mu*eta_2)/(1+2*mu*eta_1) <= -tau)
        condition_5 = abs(check) > a*tau
        if (!is.na(condition_1) && condition_1){
          updateW <- check - mu*tau
        } else if (!is.na(condition_2) && condition_2){
          updateW <- check + mu*tau
          # Next is case(ii) of SCAD penalty
        } else if (!is.na(condition_3) && condition_3){
          updateW <- (check - mu*eta_2)/(1+2*mu*eta_1)
        } else if (!is.na(condition_4) && condition_4){
          updateW <- (check + mu*eta_2)/(1+2*mu*eta_1)
          # Next is finally case (iii) of SCAD penalty
        } else if (!is.na(condition_5) && condition_5){
          updateW <- check
        } else {
          updateW <- old
        }
        finalUpdate <- updateW
    } else if (penalty == "LASSO"){
        condition_1 = check > mu*tau
        condition_2 = check < -mu*tau
        if (!is.na(condition_1) && condition_1){
            updateW <- check - mu*tau
        } else if (!is.na(condition_2) && condition_2){
            updateW <- check + updateW
        } else {
            updateW <- old
        }
        finalUpdate <- updateW
    } else if (penalty == "ELASTIC-NET"){
        condition_1 <- check > mu*tau_EN
        condition_2 <- check < -mu*tau_EN
        if (!is.na(condition_1) && condition_1){
            updateW <- (check - mu*tau_EN)/(1+2*mu*tau)
        } else if (!is.na(condition_2) && condition_2){
            updateW <- (check + mu*tau_EN)/(1+2*mu*tau)
        } else {
            updateW <- old
        }
        finalUpdate <- updateW
    } else {
        stop("Please provide the penalty function to be either LASSO, ELASTIC-NET or SCAD")
    }
    ## Check for NA
    if (any(is.na(finalUpdate))){
        finalUpdate[which(is.na(finalUpdate))] = 0
    }
    return(finalUpdate)
}

## Function for the z-update
updateZ <- function(x, old){
  check <- sqrt(sum(x^2))
  condition <- (check <= 1)
  if (!is.na(condition)){
      if (condition){
        update_z <- x
        } else {
          update_z <- x/check
      }
  } else {
      update_z = old
  }
  return(update_z)
}

relPMDCCA.algo <- function(X_1, X_2, lambda, tauW_1, tauW_2, initW_1 = NULL, initW_2 = NULL, a = 3.7, sd.mu = sqrt(2), nIter = 1000, penalty = "LASSO", element_wise = TRUE, tau_EN = NULL){
    # We use X notation for X_1 and Y for X_2
    X <- scale(X_1)
    Y <- scale(X_2)
    # mu must satisfy: 0 < mu <= lambda/norm(X)^2
    ## Each data-set has its own mu (and lambda) values - we fix lambda and have separate mu
    mu.x <- rnorm(1,mean = 0, sd = sd.mu)
    if ((mu.x < 0) | (mu.x > lambda/(Matrix::norm(X, "f")^2))){
        mu.x <- runif(1, min = 0, max = lambda/(Matrix::norm(X, "f")^2))
    }
    mu.y <- rnorm(1,mean = 0, sd = sd.mu)
    if ((mu.y < 0) | (mu.y > lambda/(Matrix::norm(Y, "f")^2))){
        mu.y <- runif(1, min = 0, max = lambda/(Matrix::norm(Y, "f")^2))
    }
    # Initialisations
    if (is.null(initW_1)){
        newW_1 = rep(0.1,ncol(X))
    } else {
        newW_1 <- initW_1
    }
    if (is.null(initW_2)){
        newW_2 = rep(0.1,ncol(Y))
    } else {
        newW_2 <- initW_2
    }
    newZ <- rep(0.4, nrow(X))
    newXi <- rep(1,nrow(X))

    # Run the algorithm
    for (i in 1:nIter){
        oldW_1 <- newW_1
        oldW_2 <- newW_2
        oldZ <- newZ
        oldXi <- newXi
        # Update First dataset
        # Compute omega:
        omega <- oldW_1 - ( (mu.x/lambda)*(t(X)%*%(X%*%oldW_1 - oldZ + oldXi)) )
        # Compute c:
        c = t(X)%*%Y%*%oldW_2
        for (j in 1:length(newW_1)){
            newW_1[j] <- updateW.relPMDCCA(omega = omega[j], mu = mu.x, c = c[j], tau = tauW_1, a = a, old = oldW_1[i], penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        }
        check_Z <- X%*%newW_1 + oldXi
        newZ <- updateZ(check_Z, oldZ)
        newXi <- oldXi + X%*%newW_1 - newZ
        # Update Second dataset
        oldW_1 <- newW_1
        oldW_2 <- newW_2
        oldZ <- newZ
        oldXi <- newXi
        # Compute omega:
        omega <- oldW_2 - ( (mu.y/lambda)*(t(Y)%*%(Y%*%oldW_2 - oldZ + oldXi)) )
        # Compute c:
        c = t(oldW_1)%*%t(X)%*%Y
        for (j in 1:length(newW_2)){
            newW_2[j] <- updateW.relPMDCCA(omega = omega[j], mu = mu.y, c = c[j], tau = tauW_2, a = a, old = oldW_2[i], penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        }
        check_Z <- Y%*%newW_2 + oldXi
        newZ <- updateZ(check_Z, oldZ)
        newXi <- oldXi + Y%*%newW_2 - newZ
    }
    return(list("w_1" = newW_1, "w_2" = newW_2, "z" = newZ, "xi" = newXi))
}

multi.relPMDCCA.algo <- function(X, lambda, tau, a = 3.7, sd.mu = sqrt(2), nIter = 1000, penalty = "LASSO", element_wise = TRUE, tau_EN = NULL){
    # Initialisations
    length.list <- length(X)
    n.samples = nrow(X[[1]])
    for (i in 1:length.list){
        X[[i]] <- scale(X[[i]])
        if (any(is.na(X[[i]]))){
            stop("Please remove/replace any missing data before implementing sCCA.")
        }
        if (nrow(X[[i]]) != n.sample){
            stop("Please make sure that all datasets in list X have the same number of samples")
        }
    }
    # mu must satisfy: 0 < mu <= lambda/norm(X)^2
    ## Each data-set has its own mu (and lambda) values - we fix lambda and have separate mu
    mu.w <- vector("list", length = length.list)
    newW <- vector("list", length = length.list)
    for (i in 1:length.list){
      mu.w[[i]] <- rnorm(1,mean = 0, sd = sd.mu)
      if ((mu.w[[i]] < 0) | (mu.w[[i]] > lambda/(Matrix::norm(X[[i]], "f")^2))){
        mu.w[[i]] <- runif(1, min = 0, max = lambda/(Matrix::norm(X[[i]], "f")^2))
      }
      newW[[i]] <- rep(0.1,ncol(X[[i]]))
    }
    newZ <- rep(0.4, nrow(X[[1]]))
    newXi <- rep(1,nrow(X[[1]]))

    # Run the algorithm
    for (i in 1:nIter){
        for (j in 1:length.list){
            oldW <- newW
            oldZ <- newZ
            oldXi <- newXi
            # Compute omega
            omega <- oldW[[j]] - ( (mu.w[[j]]/lambda)*(t(X[[j]])%*%(X[[j]]%*%oldW[[j]] - oldZ + oldXi)) )
            # Compute c:
            Yv <- 0
            for (k in 1:length.list){
                if (k !=j){
                    Yv <- Yv + X[[k]]%*%oldW[[k]]
                }
            }
            c <- t(X[[j]])%*%
            for (k in 1:length(newW[[j]])){
                newW[[j]][k] <- updateW.relPMDCCA(omega = omega[k], mu = mu.w[[j]], c = c[k], tau = tau[[j]], a = a, old = oldW[[j]][k], penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
            }
            check_Z <- X[[j]]%*%newW[[j]] + oldXi
            newZ <- updateZ(check_Z, oldZ)
            newXi <- oldXi + X[[j]]%*%newW[[j]] - newZ
        }
    }
    return(list("w" = newW, "z" = newZ, "xi" = newXi))
}

relPMDCCA <- function(X_1, X_2, lambda, tauW_1, tauW_2, initW_1 = NULL, initW_2 = NULL, a = 3.7, sd.mu = sqrt(2), nIter = 1000, penalty = "LASSO", element_wise = TRUE, tau_EN = NULL, R = 1){
    if (R == 1){
        values = relPMDCCA.algo(X_1 = X_1, X_2 = X_2, lambda = lambda, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        finalW_1 <- values$w_1
        finalW_2 <- values$w_2
    } else if (R < 1){
        stop("Please provide a positive integer for the number of canonical pairs (R)")
    } else {
        values = relPMDCCA.algo(X_1 = X_1, X_2 = X_2, lambda = lambda, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        tempW_1 <- values$w_1
        tempW_2 <- values$w_2
        for (i in 1:(R-1)){
            addX_1_1 <- t(tempW_1) %*% t(X_1) %*% X_1
            addX_1_2 <- t(tempW_2) %*% t(X_2) %*% X_1
            addX_2_1 <- t(tempW_2) %*% t(X_1) %*% X_2
            addX_2_2 <- t(tempW_1) %*% t(X_1) %*% X_2
            Xtile_1 <- rbind(X_1, addX_1_1, addX_1_2)
            Xtilde_2 <- rbind(X_2, addX_2_1, addX_2_2)
            values = relPMDCCA.algo((X_1 = Xtilde_1, X_2 = Xtilde_2, lambda = lambda, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN))
            tempW_1 <- cbind(tempW_1, values$w_1)
            tempW_2 <- cbind(tempW_2, values$w_2)
        }
        finalW_1 <- tempW_1
        finalW_2 <- tempW_2
    }
    return(list("W_1" = finalW_1, "W_2" = finalW_2))
}

multi.relPMDCCA <- function(X, lambda, tau, a = 3.7, sd.mu = sqrt(2), nIter = 1000, penalty = "LASSO", element_wise = TRUE, tau_EN = NULL, R = 1){
    if (R == 1){
        values = multi.relPMDCCA.algo(X = X, lambda = lambda, tau = tau, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        finalW <- values$w
    } else if (R < 1){
        stop("Please provide a positive integer for the number of canonical pairs (R)")
    } else {
        values = multi.relPMDCCA.algo(X = X, lambda = lambda, tau = tau, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        tempW <- values$w
        for (k in 1:(R-1)){
            Xtilde <- vector("list", length = length(X))
            addX <- vector("list", length = length(X))
            for (i in 1:length(X)){
                for (j in 1:length(X)){
                    addX[[j]] <- t(tempW[[j]])%*%t(X[[j]]) %*% X[[i]]
                }
                Xtilde[[i]] <- X[[i]]
                for (j in 1:length(X)){
                    Xtilde[[i]] <- rbind(Xtilde[[i]], addX[[j]])
                }
            }
            values = multi.relPMDCCA.algo(X = Xtilde, lambda = lambda, tau = tau, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
            for (j in 1:length(tempW)){
                tempW[[j]] <- cbind(tempW[[j]], values$w[[j]])
            }
        }
        finalW <- tempW
    }
    return(list("W" = finalW))
}
