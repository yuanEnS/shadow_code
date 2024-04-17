
# 1. P_sz -----------------------------------------------------------------
expit <- function(z){
  return (1/(1 + exp(-z)))
}
Neg_MLE_P_s1 <- function(data,par){
  ## data: c(C,A); par: initial value of xi_P_11
  data <- data[data$Z==1,]
  data$s1 <- expit(as.matrix(data[,1:(ncol(C)+1)]) %*% par)  
  data$s0 <- 1 - expit(as.matrix(data[,1:(ncol(C)+1)]) %*% par)  
  return(
      -sum(
        c(log(data[data$S==1,]$s1),log(data[data$S==0,]$s0))
      )
    )
}


Neg_MLE_P_s0 <- function(data,par){
  ## data: c(C,A); par: initial value of xi_P_10
  data <- data[data$Z==0,]
  data$s1 <- expit(as.matrix(data[,1:(ncol(C)+1)]) %*% par)
  data$s0 <- 1 - expit(as.matrix(data[,1:(ncol(C)+1)]) %*% par)  
  return(
    -sum(
      c(log(data[data$S==1,]$s1),log(data[data$S==0,]$s0))
    )
  )
}


Neg_MLE_P_Y1111 <- function(data,par){
  ## data: c(C,A); par: initial value of xi_y_1111
  data <- data[data$Z==1,]
  data <- data[data$S==1,]
  data <- data[data$R==1,]
  data$y1 <- expit(as.matrix(data[,1:(ncol(C)+1)]) %*% par)
  data$y0 <- 1 - data$y1
  return(
    -sum(
      c(log(data[data$Y==1,]$y1),log(data[data$Y==0,]$y0))
    )
  )
}


Neg_MLE_P_Y1101 <- function(data,par){
  ## data: c(C,A); par: initial value of xi_y_1101
  data <- data[data$Z==0,]
  data <- data[data$S==1,]
  data <- data[data$R==1,]
  data$y1 <- expit(as.matrix(data[,1:(ncol(C)+1)]) %*% par)
  data$y0 <- 1 - data$y1
  return(
    -sum(
      c(log(data[data$Y==1,]$y1),log(data[data$Y==0,]$y0))
    )
  )
}

Neg_MLR_R_111 <- function(data,par){
  ## data: c(C,A); par: initial value of xi_r_111
  data <- data[data$Z==1,]
  data <- data[data$S==1,]
  data$r1 <- expit(as.matrix(data[,1:(ncol(C)+1)]) %*% par)
  data$r0 <- 1 - data$r1
  return(
    -sum(
      c(log(data[data$R==1,]$r1),log(data[data$R==0,]$r0))
    )
  )
}

Neg_MLR_R_110 <- function(data,par){
  ## data: c(C,A); par: initial value of xi_r_111
  data <- data[data$Z==0,]
  data <- data[data$S==1,]
  data$r1 <- expit(as.matrix(data[,1:(ncol(C)+1)]) %*% par)
  data$r0 <- 1 - data$r1
  return(
    -sum(
      c(log(data[data$R==1,]$r1),log(data[data$R==0,]$r0))
    )
  )
}
