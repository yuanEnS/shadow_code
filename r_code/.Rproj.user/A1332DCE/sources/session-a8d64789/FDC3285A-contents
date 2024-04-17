### Parameter Estimation:
### E(Y | Z = 1, G = g, C; gamma_g) = \mu_g(C; \gamma_g) = expit(C * gamma_g)
library(BB)

# source("./1_functions/0_basic_functions.R")
gammaPred_BB <- function(ini_val_gamma_LL,
                         ini_val_gamma_LD,
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
  # The function gammaPred is used to predict parameters in E[Y | Z = 1, G = g, C] = \mu_g(C ; \gamma_g)
  
  # input:
  # ini_val_gamma_LL: initial value of gamma_LL
  # ini_val_gamma_LD: initial value of gamma_LD
  # C: observed covariates;
  # A: another observed covariate(show up in assumption)
  # Z: treatment
  # S: suvival status
  # Y: outcome
  # R: missing indicator
  # dim_status: dimensions
  # hat_alpha: estimated alpha
  # hat_beta1: estimated beta1
  # hat_beta01: estimated beta01
  
  # output:
  # the estimates of gamma_LL and gamma_LD in E[Y | Z = 1, G = g, C] = \mu_g(C ; \gamma_g)
  
  x <- cbind(C, A, Z, S, Y, R)
  
  # Moments function
  g_gamma <- function(gamma) {
    ### Moments function of gamma_LL and gamma_LD in E[Y | Z = 1, G = g, C] = \mu_g(C ; \gamma_g)
    
    # input:
    # gamma: c(gamma_LL and gamma_LD)
    # x: cbind(C,A,Z,S,Y,R)
    
    # output:
    # Moments to be equal to zero.
    len_gamma = length(gamma)
    gamma_LL <- gamma[1:(len_gamma / 2)]
    gamma_LD <- gamma[(len_gamma / 2 + 1):len_gamma]
    
    C <- x[, 1:(dim_status[1])]
    A <- x[, (dim_status[1] + 1)]
    Z <- x[, (dim_status[1] + 2)]
    S <- x[, (dim_status[1] + 3)]
    Y <- x[, (dim_status[1] + 4)]
    R <- x[, (dim_status[1] + 5)]
    
    pr_R1 <- expit(as.matrix(cbind(C, A, Y)) %*% hat_alpha)
    
    pr_G_LL_ZS1 <- expit(as.matrix(cbind(C, A)) %*% hat_beta01)
    pr_G_LD_ZS1 <- 1 - pr_G_LL_ZS1
    
    mu_LL <- expit(as.matrix(C) %*% gamma_LL)
    mu_LD <- expit(as.matrix(C) %*% gamma_LD)
    m_sum <- mu_LL * pr_G_LL_ZS1 + mu_LD * pr_G_LD_ZS1
    
    m1 <- (R / pr_R1) * (Y - m_sum) * Z * S
    
    pi_0 <- pr_G_LL_ZS1
    pi_0_c <- pi_0 * C[, 2]
    pi_1 <- pr_G_LD_ZS1
    pi_1_c <- pi_1 * C[, 2]
    
    h_dim_expand <- as.matrix(cbind(pi_0, pi_1, pi_0_c, pi_1_c))
    moments <- NULL
    for (i in 1:ncol(h_dim_expand)) {
      moments <- cbind(moments, m1 * h_dim_expand[, i])
    }
    return(c(sum(moments[,1]),sum(moments[,2]),sum(moments[,3]),sum(moments[,4])))
  }
  
  ini_gamma <- c(ini_val_gamma_LL, ini_val_gamma_LD)
  res <- BBsolve(par  = ini_gamma,fn = g_gamma,method = 3)
  return(res)
}