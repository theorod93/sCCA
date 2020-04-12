#############################################################################################################################
######################### Evaluation Measures ###############################################################################
#############################################################################################################################


## Non-zero(NZ):
computeNZ <- function(uHat){
  uHat <- round(uHat, 3)
  outTemp <- length(which(uHat != 0))
  return(outTemp)
}
## True Non-zero(TRUENZ)
computeTRUENZ <- function(uTrue, uHat){
  outTemp <- 0
  uHat <- round(uHat, 3)
  for(i in 1:length(uHat)){
    if (uTrue[i] != 0 && uHat[i] != 0){
      outTemp = outTemp + 1
    }
  }
  return(outTemp)
}
## True zero(TRUEZ)
computeTRUEZ <- function(uTrue, uHat){
  outTemp <- 0
  uHat <- round(uHat, 3)
  for(i in 1:length(uHat)){
    if (uHat[i] == uTrue[i] && uHat[i] == 0){
      outTemp = outTemp + 1
    }
  }
  return(outTemp)
}
## Angle distance(ANGLE)
computeANGLE <- function(uTrue, uHat){
  outTemp <- sqrt(abs(1-(as.numeric(t(uTrue)%*%uHat)^2)))
  return(outTemp)
}
## Loss function(LOSS)
computeLOSS <- function(uTrue, uHat){
  minVal1 <- computeNorm(uTrue-uHat)^2
  minVal2 <- computeNorm(uTrue+uHat)^2
  outTemp <- min(minVal1, minVal2)
  return(outTemp)
}
## Get norm of a vector
computeNorm <- function(x) sqrt(sum(x^2))

## Confusion Matrix for zero Vs non-zero elements
computeConfusionMatrix <- function(uTrue, uHat){
  if (length(uTrue) != length(uHat)){
    stop("Please provide two vectors of the same length(size)")
  }
  conf <- matrix(0,nrow = 2, ncol = 2)
  colnames(conf) <- c("nzTrue", "zTrue")
  rownames(conf) <- c("nzHat", "zHat")
  for (i in 1:length(uTrue)){
    if (uTrue[i] != 0){
      if (uHat[i] !=0){
        conf[1,1] <- conf[1,1] + 1
      } else {
        conf[2,1] <- conf[2,1] + 1
      }
    } else {
      if (uHat[i] != 0){
        conf[1,2] <- conf[1,2] + 1
      }else {
        conf[2,2] <- conf[2,2] + 1
      }
    }
  }
  return(conf)
}

## Various measure, such as, tpr, tnr, precision and negative preditive value
computeMeasure <- function(ConfMat){
  tnr <- ConfMat[2,2]/(ConfMat[2,2] + ConfMat[1,2])
  tpr <- ConfMat[1,1]/(ConfMat[1,1] + ConfMat[2,1])
  ppv <- ConfMat[1,1]/(ConfMat[1,1] + ConfMat[1,2])
  npv <- ConfMat[2,2]/(ConfMat[2,2] + ConfMat[2,1])
  return(list("tnr" = tnr, "tpr" = tpr, "ppv" = ppv, "npv" = npv))
}
