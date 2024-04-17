### Parameter Estimation: ### pr(R = 1| S = 1 ,A, C, Y; alpha)
# source("./1_functions/0_basic_functions.R")
alphaPred <- function(ini_val_alpha, C, A, Z, S, Y, R, dim_status) {
  # The function alphaPred is used to predict parameters in pr(R = 1| S = 1 ,A, C, Y; alpha)
  
  # input:
  # alpha: initial value;
  # C: observed covariates;
  # A: another observed covariate(show up in assumption)
  # Z: treatment
  # S: survival status
  # Y: outcome
  # R: missing indicator
  
  # output:
  # the estimates of alpha in pr(R = 1| S = 1 ,A, C, Y; alpha)
  
  x <- cbind(C, A, Z, S, Y, R)  ## better to be more data-adaptive
  ## Moments function
  
  g <- function(alpha, x) {
    C <- x[, 1:dim_status[1]]
    A <- x[, dim_status[1] + 1]
    Z <- x[, dim_status[1] + 2]
    S <- x[, dim_status[1] + 3]
    Y <- x[, dim_status[1] + 4]
    R <- x[, dim_status[1] + 5]
    pr_R1 <- expit(as.matrix(cbind(C, A, Y)) %*% alpha)
    
    h_dim_expand <- cbind(C, A, Y)
    
    m <- NULL
    for (i in 1:ncol(h_dim_expand)) {
      m <- cbind(m, as.matrix((R / pr_R1 - 1) * S * h_dim_expand[, i]))
    }
    return(m)
  }
  
  res <- gmm(g, x, ini_val_alpha, vcov  = "iid")
  return(res)
}