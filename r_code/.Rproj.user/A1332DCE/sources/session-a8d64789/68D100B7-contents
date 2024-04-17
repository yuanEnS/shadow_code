### Parameter Estimation: 
## parameters: pr( ZS = 1 | A,C ;beta1) = expit((C,A)*beta1) = \theta_1(A, C; \beta1)
## parameters: pr( G = LL | ZS = 1 A,C ) / pr(ZS = 1 | A,C) = expit((C,A)*beta01) = \theta_{0/1}(A, C; \beta2)

## Method: MLE
# source("./1_functions/0_basic_functions.R")
Neg_MLE <- function(data, par) {
  # The function Neg_MLE is to compute the -MLE and then to estimate the parameter
  
  # input:
  # data = cbind(C,A,Z,S), where, C are observed covariates; A is another observed covariate; S: survival status and Z = 1
  # par is initial parameter.
  
  # output:
  # - MLE
  len_par <- length(par)
  
  beta1 <- par[1:(len_par / 2)]
  beta01 <- par[(len_par / 2 + 1):len_par]
  
  data$theta1_value <- expit(as.matrix(data[, 1:(ncol(C) + 1)]) %*% beta1)
  
  data$theta01_value <- expit(as.matrix(data[, 1:(ncol(C) + 1)]) %*% beta01)
  
  df_sur <- data[data$S == 1, ]
  df_over <- data[data$S == 0, ]
  
  ## df_zs symbols
  df_11 <- df_sur[df_sur$Z == 1, ]
  df_01 <- df_sur[df_sur$Z == 0, ]
  df_10 <- df_over[df_over$Z == 1, ]
  df_00 <- df_over[df_over$Z == 0, ]
  
  ## pr_zs symbol
  
  # pr_set_11 <- as.matrix(df_11[,1:(ncol(df_11)-2)])%*%beta1
  pr_11 <- df_11$theta1_value
  
  # pr_set_01 <- as.matrix(df_01[,1:(ncol(df_01)-2)])%*%beta2
  pr_01 <- df_01$theta1_value * df_01$theta01_value
  # pr_set_10 <- as.matrix(df_10[,1:(ncol(df_10)-2)])%*%beta1
  pr_10 <- 1 - df_10$theta1_value
  # pr_set_00 <- as.matrix(df_00[,1:(ncol(df_00)-2)])%*%beta2
  pr_00 <- 1 - df_00$theta1_value * df_00$theta01_value
  pr_s <- c(log(pr_11), log(pr_01), log(pr_10), log(pr_00))
  res <- sum(-pr_s)
  
  return(res)
}
