### Parameter Estimation:
### E(Y | Z = 0, G = LL, C; eta) = m(C; \eta) = expit(C * eta)

# source("./1_functions/0_basic_functions.R")
# Moments function
g_eta <- function(eta, x) {
  ### Moments function of eta in E[Y | Z = 0, S = 1 , C] = m(C; \eta)
  
  # input:
  # eta: initial value of eta
  # x: cbind(C,A,Z,S,Y,R)
  
  # output:
  # Moments to be equal to zero.
  
  C <- as.data.frame(x[, 1:(dim_status[1])])
  A <- x[, (dim_status[1] + 1)]
  Z <- x[, (dim_status[1] + 2)]
  S <- x[, (dim_status[1] + 3)]
  Y <- x[, (dim_status[1] + 4)]
  R <- x[, (dim_status[1] + 5)]
  
  pr_R1 <- expit(as.matrix(cbind(C, A, Y)) %*% hat_alpha)
  
  m_Ceta <- expit(as.matrix(C) %*% eta)
  m <- (R / pr_R1) * (Y - m_Ceta) * (1 - Z) * S
  
  pr_G_LL_ZS1 <- expit(as.matrix(cbind(C, A)) %*% hat_beta01)
  pr_G_LD_ZS1 <- 1 - pr_G_LL_ZS1
  
  h_dim_expand <-
    cbind(pr_G_LL_ZS1, pr_G_LL_ZS1 * C[, 2:ncol(C)]) ### WHY pi_0 and pi_0 * C here?  NONONO, it makes no sense.
  
  # # try another choice
  # h_dim_expand <- C
  
  moments <- NULL
  for (i in 1:ncol(h_dim_expand)) {
    to_add <- m * h_dim_expand[, i]
    moments <- as.matrix(cbind(moments, to_add))
  }
  return(moments)
}
etaPred <-
  function(ini_val_eta,
           C,
           A,
           Z,
           S,
           Y,
           R,
           dim_status,
           hat_alpha,
           hat_beta1,
           hat_beta01) {
    # The function gammaPred is used to predict parameters in E[Y | Z = 0, S = 1 , C] = m(C; \eta)
    
    # input:
    # eta: initial value of eta
    # C: observed covariates;
    # A: another observed covariate(show up in assumption)
    # Z: treatment
    # S: survival status
    # Y: outcome
    # R: missing indicator
    # hat_alpha: estimated alpha
    # hat_beta1: estimated beta1
    # hat_beta01: estimated beta01
    
    # output:
    # the estimates of eta in E[Y | Z = 0, S = 1 , C] = m(C; \eta)
    
    x <- cbind(C, A, Z, S, Y, R)
    

    res <- gmm(g_eta, x, ini_val_eta, vcov = 'iid')
    return(res)
  }
