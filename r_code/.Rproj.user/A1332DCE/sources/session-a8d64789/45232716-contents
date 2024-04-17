Estimation_boot <- function(seeds,dat_obs,dim_status){
  seed_yuan = seeds

  dat_bootstrap_samples <- dat_obs[sample.int(nrow(dat_obs),replace = TRUE),]
  
  C <- dat_bootstrap_samples[,1:dim_status[1]]
  A <- dat_bootstrap_samples[,dim_status[1]+1]
  Z <- dat_bootstrap_samples[,dim_status[1]+2]
  S <- dat_bootstrap_samples[,dim_status[1]+3]
  Y <- dat_bootstrap_samples[,dim_status[1]+4]
  R <- dat_bootstrap_samples[,dim_status[1]+5]
  
  ### pr(R = 1 | S = 1 ,A, C, Y; alpha)
  ### pr(S = 1 | Z = 1, A ,C; beta1)
  ### pr(S = 1 | Z = 0, A ,C; beta01)
  ### E[Y | Z = 1, G = g, C ] = \mu_g(C; \gamma_g)
  ### E[Y | Z = 0, S = 1, C] = m(C; \eta)
  ## initial values
  
  # 
  # ini_val_alpha <- rep(0, (ncol(C)+2))
  # ini_val_beta1 <- rep(0, ncol(C) + 1)
  # ini_val_beta01 <- rep(0, ncol(C) + 1)
  # ini_val_gamma_LL <- rep(0,ncol(C))
  # ini_val_gamma_LD <- rep(0,ncol(C))
  # ini_val_eta <- rep(0,ncol(C))
  set.seed(seed_yuan)
  ini_val_alpha <- runif((ncol(C)+2))
  ini_val_beta1 <- runif(ncol(C) + 1)
  ini_val_beta01 <- runif(ncol(C) + 1)
  ini_val_gamma_LL <- runif(ncol(C))
  ini_val_gamma_LD <- runif(ncol(C))
  ini_val_eta <- runif(ncol(C))
  
  hat_alpha <-
    alphaPred(ini_val_alpha = ini_val_alpha,
              C,
              A,
              Z,
              S,
              Y,
              R,
              dim_status = dim_status)$coefficients
  
  beta_estimate <-
    optim(
      c(ini_val_beta1, ini_val_beta01),
      Neg_MLE,
      data = select(dat_bootstrap_samples, c('bias',covariates,'A','Z','S')),
      method = "BFGS"
    )$par
  hat_beta1 <- beta_estimate[1:(ncol(C) + 1)]
  hat_beta01 <- beta_estimate[(ncol(C) + 2):length(beta_estimate)]
  rm(list = c("beta_estimate"))
  
  
  hat_gamma_BB<- gammaPred_BB(ini_val_gamma_LL,
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
                              hat_beta01)
  
  hat_gamma <- hat_gamma_BB$par
  hat_gamma_LL <- hat_gamma[1:(length(hat_gamma)/2)]
  hat_gamma_LD <- hat_gamma[(length(hat_gamma)/2+1):length(hat_gamma)]
  
  ### singular:: tA^T) = rank(A), A is gradients of m_sum
  hat_eta <- etaPred(
    ini_val_eta = ini_val_eta,
    C,
    A,
    Z,
    S,
    Y,
    R,
    dim_status = dim_status,
    hat_alpha = hat_alpha,
    hat_beta1 = hat_beta1,
    hat_beta01 = hat_beta01
  )$coefficients
  
  # 3. Estimation of causal effect ------------------------------------------
  ### step 1: Estimate E(Y_0 | G = LL)
  
  # f(G = LL | A, C)
  theta1_value <-  expit(as.matrix(cbind(C, A)) %*% hat_beta1)
  theta01_value <- expit(as.matrix(cbind(C, A)) %*% hat_beta01)
  fLLAC <- theta1_value * theta01_value
  rm(list = c("theta01_value","theta1_value"))
  
  # E f(G = LL | A,C)
  ELLAC <- mean(fLLAC)
  
  ## E(Y_0 | G = LL ,A ,C) and E(Y_0 | G = LL)
  E0LLAC <- expit(as.matrix(C) %*% hat_eta)
  E0LL <- mean(E0LLAC * fLLAC) / ELLAC
  print(E0LL)
  
  ### step 2: Estimate E(Y_1 | G = LL),
  E1LLAC <- expit(as.matrix(C) %*% hat_gamma_LL)
  E1LL <- mean(E1LLAC * fLLAC) / ELLAC
  print(E1LL)
  
  ## step 3: Causal effect
  CE_cur <- E1LL - E0LL
  return(CE_cur)
}

