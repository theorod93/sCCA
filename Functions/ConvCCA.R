###############################################################################
########### ConvCCA ###########
###############################################################################

# To raise a matrix x to a given n power
"%^%" <- function(x, n){
    with(eigen(x), vectors %*% (values^n * t(vectors)))
}

## Functions to update canonical vectors within convCCA algorithm
# penalty can take "SCAD", or "LASSO" all capital
updateW.convCCA <- function(K, v, tau, penalty = "LASSO"){
    # 4-step process
    # Step 1 -- Initialise
    K <- as.matrix(K)
    v <- as.numeric(v)
    vx <- K %*% v

    # Step 2 -- Normalise
    norm.vx <- as.numeric(sqrt(t(vx)%*%vx))
    if(norm.vx==0) norm.vx <- 1
    vx <- vx / norm.vx

    # Step 3 -- Implement algorithm via penalty function
    names(vx) <- 1:length(vx)
    if (penalty == "SCAD"){
        # First part
        vx1 <- vx[abs(vx)<=2*tau]
        u.new1 <- abs(vx1) - tau
        u.new1 <- (u.new1 + abs(u.new1))/2
        u.new1 <- u.new1*sign(vx1)
        # Second part
        vx2 <- vx[abs(vx)>2*tau & abs(vx)<=3.7*tau]
        u.new2 <- (2.7*vx2-sign(vx2)*3.7*tau)/1.7
        # Third part
        vx3 <- vx[abs(vx)>3.7*tau]
        u.new3 <- vx3
        u.new <- c(u.new1,u.new2,u.new3) # Combine
        u.new <- u.new[names(vx)] # Puts the data back to the original order
        u.new <- as.numeric(u.new) # Remove the names
    } else if (penalty == "LASSO"){
        u.new <- abs(vx) - tau
        u.new <- (u.new + abs(u.new))/2
        u.new <- u.new*sign(vx)
        u.new <- as.numeric(u.new)
    } else {
        stop("Please provide the penalty function to be either LASSO or SCAD")
    }

    # Step 4 -- Normalise again
    norm.u.new <- as.numeric(sqrt(t(u.new)%*%u.new))
    if(norm.u.new==0) norm.u.new <- 1
    u.new <- u.new / norm.u.new
    return(u.new)
}

updateWmulti.convCCA <- function(K, w, whichDset, tau, penalty = "LASSO"){
    if (class(K) != "list"){
        stop("Please provide a valid K Correlation matrix. It should be a list")
    }

    # 4-step process
    # Step 1 -- Initialise
    vx <- 0
    # Take the sum of canonical pairs to use in updating w
    for (i in 1:length(K)){
        if (i != whichDSet){
            K[[i]] <- as.matrix(K[[i]])
            w[[i]] <- as.numeric(w[[i]])
            vx <- vx + K[[i]] %*% w[[i]]
        }
    }

    # Step 2 -- Normalise
    norm.vx <- as.numeric(sqrt(t(vx)%*%vx))
    if (is.na(norm.vx)) norm.vx <- 1
    if(norm.vx==0) norm.vx <- 1
    vx <- vx / norm.vx

    # Step 3 - Implement penalty Functions
    if (penalty == "SCAD"){
        # First part
        vx1 <- vx[abs(vx)<=2*tau]
        u.new1 <- abs(vx1) - tau
        u.new1 <- (u.new1 + abs(u.new1))/2
        u.new1 <- u.new1*sign(vx1)
        # Second part
        vx2 <- vx[abs(vx)>2*tau & abs(vx)<=3.7*tau]
        u.new2 <- (2.7*vx2-sign(vx2)*3.7*tau)/1.7
        # Third part
        vx3 <- vx[abs(vx)>3.7*tau]
        u.new3 <- vx3
        u.new <-  c(u.new1,u.new2,u.new3)
        u.new <- u.new[names(vx)] # Puts the data back to the original order
        u.new <- as.numeric(u.new) # Remove the names
    } else if (penalty == "LASSO"){
        u.new <- abs(vx) - tau
        u.new <- (u.new + abs(u.new))/2
        u.new <- u.new*sign(vx)
        u.new <- as.numeric(u.new)
    } else {
        stop("Please provide the penalty function to be either LASSO or SCAD")
    }

    # Step 4 -- Normalise again
    norm.u.new <- as.numeric(sqrt(t(u.new)%*%u.new))
    if (is.na(norm.u.new)) norm.u.new <- 1
    if(norm.u.new==0) norm.u.new <- 1
    u.new <- u.new / norm.u.new
    return(u.new)
}

gettingInverse <- function(X){
    A <- X
    svdA <- svd(A)
    d <- round(svdA$d,3)
    or <- order(d, decreasing = TRUE)
    U <- svdA$u[or,]
    V <- svdA$v[,or]
    d <- d[or]
    l <- length(which(d == 0))
    val <- nrow(A)-l
    S <- diag(svdA$d)
    Splus <- diag(c(1/svdA$d[1:val],rep(0,l)))
    Aplus <- V%*%Splus%*%t(U)
    return(Aplus)
}

# Main function for ConvCCA
convCCA.algo <- function(X_1, X_2, tauW_1, tauW_2, initW_1 = NULL, initW_2 = NULL, nIter = NULL, penalty = "LASSO"){
    # Initialisations
    # We use X notation for X_1 and Y for X_2
    X <- scale(X_1)
    Y <- scale(X_2)
    SigmaX <- cov(X)
    SigmaY <- cov(Y)
    SigmaXY <- cov(X, Y)
    if (ncol(X) > nrow(X)){
        SigmaXforK <- diag(ncol(X))
    } else {
        SigmaXforK <- SigmaX %^% (-0.5)
    }
    if (ncol(Y) > nrow(Y)){
        SigmaYforK <- diag(ncol(Y))
    } else {
        SigmaYforK <- SigmaY %^% (-0.5)
    }
    # Compute K
    K <- SigmaXforK %*% SigmaXY %*% SigmaYforK
    if (!any(is.na(K))){
        ee <- eigen(t(K)%*%K)
        d <- ee$values[1]
    }

    # Initialise
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

    # Run the algorithm:
    if (is.null(nIter)){
        diff <- 1
        diffTemp <- 1
        while (diff >1e-05 && !is.na(diffTemp)) {
            oldW_1 <- newW_1
            oldW_2 <- newW_2
            ## Update u:
            newW_1 <- updateW.convCCA(K = K, v = oldW_2, tau = tauW_1, penalty = penalty)
            ## Update v:
            newW_2 <- updateW.convCCA(K = t(K), v = newW_1, tau = tauW_2, penalty = penalty)
            diff0.v <- sqrt(sum((oldW_2 - newW_2) ^ 2))
            diff0.u <- sqrt(sum((oldW_1 - newW_1) ^ 2))
            diffTemp <- diff0.u + diff0.v
            if (!is.na(diffTemp)){
                diffU = diffTemp
            }
        }
    } else if (!is.numeric(nIter)){
        stop("Please provide a numerical value for the number of iterations(nIter)")
    } else {
        for (i in 1:nIter){
          oldW_1 <- newW_1
          oldW_2 <- newW_2
          newW_1 <- updateW.convCCA(K = K, v = oldW_2, tau = tauW_1, penalty = penalty)
          newW_2 <- updateW.convCCA(K = t(K), v = newW_1, tau = tauW_2, penalty = penalty)
        }
    }
    return(list("w_1" = newW_1, "w_2" = newW_2, "d" = d))
}

# Main function for multi ConvCCA
multi.convCCA.algo <- function(X, tau, nIter = 1000, penalty = "LASSO", initW = NULL){
    # Initialisations
    length.list <- length(X)
    s <- nrow(X[[1]])
    for (i in 1:length.list){
        X[[i]] <- scale(X[[i]])
        if (any(is.na(X[[i]]))){
            stop("Please remove/replace any missing data before implementing CCA.")
        }
        if (nrow(X[[i]]) != n.samples){
            stop("Please make sure that all data-sets in X have the same number of samples.")
        }
    }
    wNew <- vector("list", length = length.list)
    wOld <- vector("list", length = length.list)
    SigmaX <- vector("list", length = length.list)
    SigmaXforK <- vector("list", length = length.list)
    if (is.null(initW)){
        for (i in 1:length.list){
            wNew[[i]] <- rep(0.1,ncol(X[[i]]))
            wOld[[i]] <- wNew[[i]]
            SigmaX[[i]] <- cov(X[[i]])
            if (ncol(X[[i]]) > nrow(X[[i]])){
                SigmaXforK[[i]] <- diag(ncol(X[[i]]))
            } else {
                SigmaXforK[[i]] <- SigmaX[[i]] %^% (-0.5)
            }
        }
    } else {
        for (i in 1:length.list){
          wNew[[i]] <- initW[[i]]
          wOld[[i]] <- wNew[[i]]
          SigmaX[[i]] <- cov(X[[i]])
          if (ncol(X[[i]]) > nrow(X[[i]])){
              SigmaXforK[[i]] <- diag(ncol(X[[i]]))
          } else {
              SigmaXforK[[i]] <- SigmaX[[i]] %^% (-0.5)
          }
        }
    }
    diff <- 1
    n.t <- 0
    # Run the algorithm
    for (k in 1:nIter){
        for (i in 1:length.list){
            # Compute K matrices
            K <- vector("list", length = length.list-1)
            for (j in 1:length.list){
                if (j!=i){
                    K[[j]] = SigmaXforK[[i]] %*% cov(X[[i]], X[[j]])%*%SigmaXforK[[j]]
                }
            }
            wOld[[i]] <- wNew[[i]]
            ## Update u:
            wNew[[i]] <- updateWmulti.convCCA(K = K, w = wOld, whichDSet = i, tau = tau[[i]], penalty = penalty)
        }
    }
    return(list("w" = wNew))
}

convCCA <-  function(X_1, X_2, tauW_1, tauW_2, initW_1 = NULL, initW_2 = NULL, nIter = NULL, penalty = "LASSO", R = 1){
    if (R == 1){
        values = convCCA.algo(X_1 = X_1, X_2 = X_2, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, nIter = nIter, penalty = penalty)
        finalW_1 <- values$w_1
        finalW_2 <- values$w_2
    } else if (R < 1){
        stop("Please provide a positive integer for the number of canonical pairs (R)")
    } else {
        values = convCCA.algo(X_1 = X_1, X_2 = X_2, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, nIter = nIter, penalty = penalty)
        tempW_1 <- values$w_1
        tempW_2 <- values$w_2
        for (i in 1:R-1){
            addX_1_1 <- t(tempW_1) %*% t(X_1) %*% X_1
            addX_1_2 <- t(tempW_2) %*% t(X_2) %*% X_1
            addX_2_1 <- t(tempW_2) %*% t(X_1) %*% X_2
            addX_2_2 <- t(tempW_1) %*% t(X_1) %*% X_2
            Xtile_1 <- rbind(X_1, addX_1_1, addX_1_2)
            Xtilde_2 <- rbind(X_2, addX_2_1, addX_2_2)
            values = convCCA.algo(Xtilde_1 = Xtilde_1, X_2 = X_2, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, nIter = nIter, penalty = penalty)
            tempW_1 <- values$w_1
            tempW_2 <- values$w_2
        }
        finalW_1 <- tempW_1
        finalW_2 <- tempW_2
    }
    return(list("W_1" = finalW_1, "W_2" = finalW_2))
}

multi.convCCA <- function(X, tau, nIter = 1000, penalty = "LASSO", initW = NULL, R = 1){
    if (R == 1){
        values = multi.convCCA.algo(X = X, tau = tau, nIter = nIter, penalty = penalty, initW = initW)
        finalW <- values$w
    } else if (R < 1){
        stop("Please provide a positive integer for the number of canonical pairs (R)")
    } else {
        values = multi.convCCA.algo(X = X, tau = tau, nIter = nIter, penalty = penalty, initW = initW)
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
            values = multi.convCCA.algo(X = Xtilde, tau = tau, nIter = nIter, penalty = penalty, initW = initW)
            for (j in 1:length(tempW)){
                tempW[[j]] <- cbind(tempW[[j]], values$w[[j]])
            }
        }
        finalW <- tempW
    }
    return(list("W" = finalW))
}
