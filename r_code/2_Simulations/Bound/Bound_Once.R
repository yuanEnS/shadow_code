library(MASS)
library(stats)
library(dplyr)
library(gmm)  
library(BB)
library(Rlab)
true_causal <- 0.3322807


# 1. generating data ------------------------------------------------------
seeds <- 20230407
seeds <- 20240426
seeds <- 20240507
data_generate <- function(seeds)
{
  
  expit <- function(linear_term)
  {
    ### expit function
    ## input: linear term
    ## output: 1/(1+exp(-linear_term))
    
    return(1/(1+exp(-linear_term)))
  }
  # 1. generating data ------------------------------------------------------
  
  set.seed(seeds)
  N <- 20000
  
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
  para_pos_C <- 0.3
  para_pos_A <- 0.7
  
  
  ## covariates setting
  C1 <- rnorm(N, mean = para_mean_C, sd = para_sd_C)

  C <- cbind(rep(1, N), C1)
  A <- rnorm(N, mean = para_mean_A, sd = para_sd_A)
  
  
  Z <- rbinom(N, 1, para_mean_Z)
  ## princial strata G
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
    Y_1_ll[i] <- rbinom(1, 1, p_LL[i])
    Y_0_ll[i] <- rbinom(1, 1, p_eta[i])
    Y_1_ld[i] <- rbinom(1, 1, p_LD[i])
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


expit <- function(linear_term)
{
  ### expit function
  ## input: linear term
  ## output: 1/(1+exp(-linear_term))
  
  return(1/(1+exp(-linear_term)))
}



dat<- data_generate(seeds)

dat$C1 <- ifelse(dat$C1>0,1,0)
dat$A <- ifelse(dat$A > 0,1,0)


# 2, Compute Some Terms ---------------------------------------------------

##### Some generic functions
f <- function(a,c,dat){
  ## compute the empirical probability of A and C
  dat_cur <- dat[dat$A == a,]
  dat_cur <- dat_cur[dat_cur$C1 == c,]
  return(dim(dat_cur)[1]/dim(dat)[1])
}
Pr_S_ACZ <- function(s,a,c,z,dat){
  ## compute the empirical probability of Pr(S = s | A = a, C = c, Z = z)
  dat_cur <- dat[(dat$A == a) & (dat$C1 == c) & (dat$Z == z),]
  dat_cur_S <- dat_cur[dat_cur$S == s, ]
  return(dim(dat_cur_S)[1]/dim(dat_cur)[1])
}
Pr_Y_SACZR <- function(y,s,a,c,z,r,dat){
  ## compute the empirical probability of Pr(Y = y | S = s, A = a,C = c,Z = z, R = r)
  dat_cur <- dat[(dat$S == s) & (dat$A == a) & (dat$C1 == c) & (dat$Z == z) & (dat$R == r),]
  dat_cur_Y <- dat_cur[dat_cur$Y == y,]
  return(dim(dat_cur_Y)[1]/dim(dat_cur)[1])
}

Pr_R_SACZ <- function(r,s,a,c,z,dat){
  dat_cur <- dat[(dat$S == s) & (dat$A == a) & (dat$C1 == c) & (dat$Z == z),]
  dat_cur_R <- dat_cur[dat_cur$R == r,]
  return(dim(dat_cur_R)[1]/dim(dat_cur)[1])
}

int_S_Z <- function(s,z,dat){
  ## compute \int_{A,C} Pr{S = s | Z = z, A, C}f(A,C) d\mu(A)d\mu(C)
  integration <- 0
  for(a in 0:1){
    for(c in 0:1){
      integration = integration + Pr_S_ACZ(s,a,c,z,dat)*f(a,c,dat)
    }
  }
  dat_cur <- dat[dat$Z ==0,]
  dat_final <- dat_cur[dat_cur$S==1,]
  res <- dim(dat_final)[1]/dim(dat_cur)[1]
  return(res)
}

int_Y_SZR <- function(y,s,z,r,dat){
  ## compute \int_{A,C} Pr(Y = y | S = s, Z = z, R =r ,A, C)f(A,C)d\mu(A)d\mu(C)
  integration <- 0
  for(a in 0:1){
    for(c in 0:1){
      integration = integration + Pr_Y_SACZR(y,s,a,c,z,r,dat)*f(a,c,dat)
    }
  }
  return(integration)
}

f_AC_LL <- function(a,c,dat){
  ## compute f(A,C | G = LL)
  res <- Pr_S_ACZ(1,a,c,0,dat)*f(a,c,dat)/int_S_Z(1,0,dat)
  return(res)
}

f_AC_S1 <- function(a,c,dat){
  ## compute f(A,C | S(1) = 1)
  res <- Pr_S_ACZ(1,a,c,1,dat)*f(a,c,dat)/int_S_Z(1,1,dat)
  return(res)
}
##### Compute values
# gamma_proportion <- int_S_Z(1,0,dat)/(1 - int_S_Z(0,1,dat))
gamma_proportion <- function(a,c,dat){
  res <- Pr_S_ACZ(1,a,c,0,dat)/Pr_S_ACZ(1,a,c,1,dat)
}
gamma_proportion2 <- function(a,c,dat){
  gamma_x <- mean(S[Z == 0 & A== a & C == c])/mean(S[Z == 1 & A== a & C == c])
  return(gamma_x)
} ## gamma_: no problem


delta_zx2 <- function(z,a,c,dat){
  res <-  mean(R[Z == z &  S == 1 & A== a & C == c])
  return(res)
}
xi_z1_x <- function(z,a,c,dat){
  res <- mean(Y[S == 1 & Z==z & R ==1 & A==a & C==c])
  return(res)
}

##### change code start
## a_u
int_a_u <- 0
for(a in 0:1){
  for(c in 0:1){
    pi_x_u <- (Pr_Y_SACZR(1,1,a,c,1,1,dat)*Pr_R_SACZ(1,1,a,c,1,dat)+Pr_R_SACZ(0,1,a,c,1,dat))
    x_gamma <- gamma_proportion(a,c,dat)
    x_a_u <- pi_x_u/x_gamma 
    int_a_u = int_a_u + min(1, x_a_u)*f_AC_LL(a,c,dat)
  }
}


## a_l
int_a_l <- 0
for(a in 0:1){
  for(c in 0:1){
    x_gamma <- gamma_proportion(a,c,dat)
    x_gamma2 <- gamma_proportion2(a,c,dat)
    pi_x_l <- Pr_Y_SACZR(1,1,a,c,1,1,dat)*Pr_R_SACZ(1,1,a,c,1,dat)
    x_a_l <- pi_x_l/x_gamma - (1 - x_gamma)/x_gamma
    to_add <- max(0,x_a_l)*f_AC_LL(a,c,dat)
    int_a_l = int_a_l + to_add 
  }
}

int_a_l <- max(0,int_a_l)



#b_u
int_b_u <- 0
for(a in 0:1){
  for(c in 0:1){
    to_add <- (Pr_Y_SACZR(1,1,a,c,0,1,dat)*Pr_R_SACZ(1,1,a,c,0,dat) + Pr_R_SACZ(0,1,a,c,0,dat))
    int_b_u <- int_b_u + to_add*f_AC_LL(a,c,dat)
  }
}

int_b_l <- 0
for(a in 0:1){
  for(c in 0:1){
    int_b_l <- int_b_l + (Pr_Y_SACZR(1,1,a,c,0,1,dat)*Pr_R_SACZ(1,1,a,c,0,dat))*f_AC_LL(a,c,dat)
  }
}

print(paste("CE falls in interval: [",int_a_l - int_b_u,int_a_u - int_b_l,"]"))

##### change code end

##### bound under randomized trial: start

# delta_z: pr(R(z) = 1 | S(z) = 1 ) = pr(R = 1 | S = 1, Z = z)
delta_z <- function(z,dat){
  dat_cur <- dat[dat$S == 1,]
  dat_cur <- dat_cur[dat_cur$Z == z,]
  dat_final <- dat_cur[dat_cur$R ==1, ]
  res <- dim(dat_final)[1]/dim(dat_cur)[1]
  return(res)
}
delta_1 <- delta_z(1,dat)

delta_0 <- delta_z(0,dat)

# Pr(S(0) = 1 | S(1) = 1) = Pr(S(0) = 1, S(1) = 1)/Pr(S(1) = 1) = Pr(S(0) = 1)/Pr(S(1) = 1) 
# =  Pr(S = 1 | Z = 0)/Pr(S = 1 | Z = 1) under randomized trial
gamma_com <- function(dat){
  dat_z_0 <- dat[dat$Z==0,]
  dat_cur_0 <- dat_z_0[dat_z_0$S==1,]
  
  dat_z_1 <- dat[dat$Z==1,]
  dat_cur_1 <- dat_z_1[dat_z_1$S==1,]
  numerator <- dim(dat_cur_0)[1]/dim(dat_z_0)[1]
  denome <- dim(dat_cur_1)[1]/dim(dat_z_1)[1]
  res <- numerator / denome
  return(res)
}

gamma <- gamma_com(dat)

# \xi_{z1} = pr{Y(z) = 1| S(z) = 1, R(z) = 1}
xi_z_1 <- function(z,dat){
  dat_cur <- dat[dat$R == 1,]
  dat_cur <- dat_cur[dat_cur$Z == z,]
  dat_cur <- dat_cur[dat_cur$S==1,]
  dat_final <- dat_cur[dat_cur$Y==1,]
  res <- dim(dat_final)[1]/dim(dat_cur)[1]
  return(res)
}
xi_01 <- xi_z_1(0,dat)
xi_11 <- xi_z_1(1,dat)

pi_1_u <- delta_1 * xi_11 + 1 - delta_1
pi_1_l <- delta_1 * xi_11 

theta_111_l <- max(0, pi_1_l/gamma - (1 - gamma)/gamma)
theta_111_u <- min(1,pi_1_u/gamma)

theta_011_u <- delta_0 * xi_01 + 1 - delta_0


theta_011_l <- delta_0 * xi_01 

print(paste("CE falls in interval: [",int_a_l - int_b_u,int_a_u - int_b_l,"]"))

print(paste("CE falls in interval: [",theta_111_l - theta_011_u,theta_111_u - theta_011_l,"]"))
##### bound under randomized trial: end
