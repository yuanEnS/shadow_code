# 1. generating data ------------------------------------------------------
data_generate <- function(seeds)
{
  
  library(MASS)
  library(stats)
  library(dplyr)
  library(gmm)
  
  
  expit <- function(linear_term)
  {
    ### expit function
    ## input: linear term
    ## output: 1/(1+exp(-linear_term))
    
    return(1/(1+exp(-linear_term)))
  }
  # 1. generating data ------------------------------------------------------
  
  set.seed(seeds)
  N <- 8000
  
  ###### some parameters ######
  
  para_delta2 <-
    0.5 ## parameters in generating beta1 and beta2 (relationships between C)
  para_beta1 <- c(1.1, para_delta2, 2) / 2
  para_beta01 <- c(0.9, -2 * para_delta2, 2) / 2
  para_mean_C <- 0
  para_sd_C <- 0.5
  para_mean_A <- 0
  para_sd_A <- 0.5
  para_mean_Z <- 0.6
  para_gamma_LL <-
    c(0.9, 0.3) # E(Y_1 | G = LL, C; gamma_LL) = expit(C * gamma_LL)
  para_gamma_LD <-
    c(0.5, 0.4) # E(Y_1 | G = LD, C; gamma_LD) = expit(C * gamma_LD)
  para_eta <- c(-0.5, 0.3) # E(Y_0 | G = LL, C; eta) = expit(C * eta)
  para_delta3 <- 0 ## parameters in missing indicator R
  para_alpha <-
    c(3, -3 * para_delta3, 1, 2.2) / 2 ## parameters in estimating indicator probability: pr(R = 1 | C,A,Y,S= 1) = expit((C,A,Y)*para_alpha)
  ## covariates setting
  C <- cbind(rep(1, N), rnorm(N, mean = para_mean_C, sd = para_sd_C))
  A <- rnorm(N, mean = para_mean_A, sd = para_sd_A)
  Z <- rbinom(N, 1, para_mean_Z)
  ## princial strata G
  para_beta1
  ## parameters: pr(ZS = 1 | A,C ;beta1) = expit((C,A)*beta1)
  para_beta01
  ## parameters: pr( G = LL|ZS = 1 A,C ) / pr(ZS = 1 | A,C) = expit((C,A)*beta01)
  
  expit_1 <-
    expit(as.matrix(cbind(C, A)) %*% para_beta1)  ## pr(ZS = 1 | A,C ;beta1)
  expit_01 <-
    expit(as.matrix(cbind(C, A)) %*% para_beta01)  ##  pr( G = LL|ZS = 1 A,C ) / pr(ZS = 1 | A,C)
  pr_G_ll <- expit_1 * expit_01
  pr_G_ld <- expit_1 * (1 - expit_01)
  pr_G_dd <- 1 - pr_G_ll - pr_G_ld
  G <- c()
  for (i in 1:N)
  {
    ## G = 0 : LL ; G = 1: LD ; G = 2 : DD
    G[i] <-
      sample(0:2,
             1,
             replace = TRUE,
             prob = c(pr_G_ll[i], pr_G_ld[i], pr_G_dd[i]))
    
  }
  rm(list = c("expit_01", "expit_1", "pr_G_ll", "pr_G_ld", "pr_G_dd"))
  
  ## survival status
  S <- rep(0, N)
  for (i in 1:N) {
    if (G[i] == 0 || (G[i] == 1 && Z[i] == 1)) {
      S[i] <- 1
    }
  }
  
  ## potential outcomes
  
  print("parameters are:")
  print(paste("gamma_LL:", para_gamma_LL[1], para_gamma_LL[2]))
  print(paste("gamma_LD:", para_gamma_LD[1], para_gamma_LD[2]))
  print(paste("eta:", para_eta[1], para_eta[2]))
  
  Y_1_ll <- c()
  Y_0_ll <- c()
  Y_1_ld <- c()
  c_multi_gamma_ll <- as.matrix(C) %*% para_gamma_LL
  c_multi_eta <- as.matrix(C) %*% para_eta
  c_multi_gamma_ld <- as.matrix(C) %*% para_gamma_LD
  p_LL <- expit(c_multi_gamma_ll)
  p_eta <- expit(c_multi_eta)
  p_LD <- expit(c_multi_gamma_ld)
  for (i in 1:N)
  {
    Y_1_ll[i] <- rbinom(1, 1, p_LL)
    Y_0_ll[i] <- rbinom(1, 1, p_eta)
    Y_1_ld[i] <- rbinom(1, 1, p_LD)
  }
  rm(
    list = c(
      "c_multi_gamma_ll",
      "c_multi_eta",
      "c_multi_gamma_ld",
      "p_LL",
      "p_eta",
      "p_LD"
    )
  )
  
  ## observational outcome Y
  Y <- c()
  for (i in 1:N)
  {
    if (G[i] == 0 && Z[i] == 1) {
      Y[i] <- Y_1_ll[i]
    }
    else if (G[i] == 0 && Z[i] == 0) {
      Y[i] <- Y_0_ll[i]
    }
    else if (G[i] == 1 && Z[i] == 1) {
      Y[i] <- Y_1_ld[i]
    } else{
      Y[i] <- -1  ## -1 to indicate survival censor
    }
  }
  
  ## missing indicator R
  R <- c()
  print("paramters used in R, alpha = ")
  print(para_alpha)
  alpha_multi_cay <- cbind(C, A, Y) %*% para_alpha
  for (i in 1:N) {
    R[i] <- rbinom(1, 1, expit(alpha_multi_cay[i]))
  }
  rm(list = c("alpha_multi_cay"))
  
  
  ## check the data generation
  dat <- cbind(C, A, Z, G, S, Y_1_ll, Y_0_ll, Y_1_ld, Y, R)
  dat <- data.frame(dat)
  names(dat) <-
    c('bias',
      'C1',
      'A',
      'Z',
      'G',
      'S',
      'Y_1_ll',
      'Y_0_ll',
      'Y_1_ld',
      'Y',
      'R')
  dim_status <-
    c(ncol(C), ncol(A), ncol(Z), ncol(S), ncol(Y), ncol(R))
  dat_obs <- select(dat, c('bias', 'C1', 'A', 'Z', 'S', 'Y', 'R'))
  
  return(dat)
}