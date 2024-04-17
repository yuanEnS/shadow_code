getwd()
## setwd("/home/u2021103709/shadow/r_code/")
library(parallel)

clnum <- detectCores()
print(paste("detectCores:" ,clnum))
clnum <- 56
print(paste("specified cores:" , clnum))
print(clnum)
cl <- makeCluster(getOption('cl.cores',clnum))
print("makecluster starts")
time_start <- Sys.time()
cyc <- 1000
ELL <- function(seeds){
  print(seeds)
  library(MASS)
  library(stats)
  library(dplyr)
  library(gmm)  
  library(BB)
  print(seeds)
  recycle_times <- 200
  true_causal <- 0.3322807

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
    # 1. generating data -------------------------------------------1.0858-----------
    
    set.seed(seeds)
    N <- 5000
    
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
    Y_1_ll <- rbinom(N, 1, p_LL)
    Y_0_ll <- rbinom(N, 1, p_eta)
    Y_1_ld <- rbinom(N, 1, p_LD)
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
    
    return(dat_obs)
  }
  
  expit <- function(linear_term)
  {
    ### expit function
    ## input: linear term
    ## output: 1/(1+exp(-linear_term))
    
    return(1/(1+exp(-linear_term)))
  }
  
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
      pi_1_c <- pi_1 * C[,2]
      
      h_dim_expand <- as.matrix(cbind(pi_0, pi_1, pi_0_c, pi_1_c))
      moments <- NULL
      for (i in 1:ncol(h_dim_expand)) {
        moments <- cbind(moments, m1 * h_dim_expand[, i])
      }
      return(c(sum(moments[,1]),sum(moments[,2]),sum(moments[,3]),sum(moments[,4])))
    }
    
    ini_gamma <- c(ini_val_gamma_LL, ini_val_gamma_LD)
    res <- BBsolve(par  = ini_gamma,fn = g_gamma)
    return(res)
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
      
      # Moments function
      g_eta <- function(eta, x) {
        ### Moments function of eta in E[Y | Z = 0, S = 1 , C] = m(C; \eta)
        
        # input:
        # eta: initial value of eta
        # x: cbind(C,A,Z,S,Y,R)
        
        # output:
        # Moments to be equal to zero.
        
        C <- x[, 1:(dim_status[1])]
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
          cbind(pr_G_LL_ZS1, pr_G_LL_ZS1 * C[, 2]) ### WHY pi_0 and pi_0 * C here?  NONONO, it makes no sense.
        
        # # try another choice
        h_dim_expand <- C
        
        moments <- NULL
        for (i in 1:ncol(h_dim_expand)) {
          to_add <- m * h_dim_expand[, i]
          moments <- as.matrix(cbind(moments, to_add))
        }
        return(moments)
      }
      res <- gmm(g_eta, x, ini_val_eta, vcov = 'iid')
      return(res)
    }
  
  
  
  dat<- data_generate(seeds)
  dat_obs <- dat[,c('bias', 'C1', 'A', 'Z', 'S', 'Y', 'R')]
  # dat_obs <- select(dat, c('bias', 'C1', 'A', 'Z', 'S', 'Y', 'R'))
  
  ######## Different from original method
  dat_obs <- dat_obs[dat_obs$R==1,]
  # dat_obs <- dat_obs[dat_obs$S==1,]
  
  C <- dat_obs[,2] # because there has existed a bias term in the package
  A <- dat_obs[,3]
  Z <- dat_obs[,4]
  S <- dat_obs[,5]
  Y <- dat_obs[,6]
  R <- dat_obs[,7]
  

# 2. estimation -----------------------------------------------------------
  library(tbd)
  sace <- sace(Z,S,Y,C,A,hessian = FALSE)
  sace_value <- sace$sace
  sace.boot <- boot.ci(sace, print.progress = FALSE)
  sace.boot$boot.sd
  sace.boot$ci
  coverage <- 1
  if(sace.boot$ci[1] > true_causal | sace.boot$ci[2] < true_causal){
    coverage <- 0
  }
# 5. return results -------------------------------------------------------
  return(c(sace_value,coverage))
}
res <- parLapply(cl,1:cyc,ELL)
CEs <- NULL
for(i in 1:cyc){
  CEs <- c(CEs,res[[i]][1])
}
coverages <- NULL
for(i in 1:cyc){
  coverages <- c(coverages,res[[i]][2])
}
time_end <- Sys.time()
time <- time_end - time_start
print(time)
save.image(paste("Repeat",cyc,"_N",5000,".RData",sep = ''))
