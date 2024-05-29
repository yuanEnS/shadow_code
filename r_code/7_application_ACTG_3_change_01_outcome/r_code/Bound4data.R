# 1. Packages -------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(gmm)
library(parallel)
library(speff2trial)

# 2. Data processing ------------------------------------------------------
data("ACTG175")
source("../1_functions/0_basic_functions.R")
source("../1_functions/1_alphaPred.R")
source("../1_functions/2_betaPred.R")
source("../1_functions/3_1_gammaPred_for_application.R")
source("../1_functions/4_etaPred_for_application.R")
source("../1_functions/5_application.R")


source("./r_code/functions4data/MLEs.R")

covariates <- c("wtkg")

median_wtkg <- median(select(ACTG175,'wtkg')$wtkg)
C <- ifelse((select(ACTG175,'wtkg'))$wtkg>median_wtkg,1,0)
# C <- cbind(C,ifelse((select(ACTG175,'karnof'))$karnof>95,1,0))
C <- as.data.frame(C)
colnames(C) <- covariates
# for(c in covariates){
#   print(c)
#   if(c == 'age' | c == 'wtkg' | c == 'karnof'){
#     C <- cbind(C,scale(select(ACTG175,c)))
#   }
#   else{
#     C <- cbind(C,select(ACTG175,c))
#   }
# }
median_A <- median(ACTG175$cd40)
A <- ifelse(ACTG175$cd40 > median_A,1,0)


Z <- (ACTG175$arms > 0)
S <- (ACTG175$cens & (ACTG175$days<7*96))
S[S] <- 1
# S <- ACTG175$offtrt==0
delta_cd4 <- ACTG175$cd496 - ACTG175$cd40
R <- (ACTG175$r == 1)
R[R] <- 1
Y <- c()
for(i in 1:length(delta_cd4)){
  if(is.na(delta_cd4[i])){
    Y <- c(Y,NA)
  }
  # else if(delta_cd4[i] - medi_delta > 0){
  else if(delta_cd4[i] > 0){
    Y <- c(Y,1)
  }
  else{
    Y <- c(Y,0)
  }
}
Y[is.na(Y)] <- -1

dat <- cbind(C,A,Z,S,Y,R)
head(dat)

# Compute Bound -----------------------------------------------------------

# 2, Compute Some Terms ---------------------------------------------------
##### Some generic functions
f <- function(a,wtkg,dat){
  ## compute the empirical probability of A and C
  dat_cur <- dat[dat$A == a,]
  dat_cur <- dat_cur[dat_cur$wtkg == wtkg,]
  res <- dim(dat_cur)[1]/dim(dat)[1]
  if(is.na(res)){
    res <- 0
    print("nan f")
  }
  return(res)
}
Pr_S_ACZ <- function(s,a,wtkg,z,dat){
  ## compute the empirical probability of Pr(S = s | A = a, C = c, Z = z)
  dat_cur <- dat[(dat$A == a) &  (dat$Z == z),]
  dat_cur <- dat_cur[dat_cur$wtkg == wtkg,]
  dat_cur_S <- dat_cur[dat_cur$S == s, ]
  res <- dim(dat_cur_S)[1]/dim(dat_cur)[1]
  if(is.na(res)){
    res <- 0
    print("NaN Pr_S_ACZ")
  }
  return(res)
}
Pr_Y_SACZR <- function(y,s,a,wtkg,z,r,dat){
  ## compute the empirical probability of Pr(Y = y | S = s, A = a,C = c,Z = z, R = r)
  dat_cur <- dat[(dat$S == s) & (dat$A == a) & (dat$Z == z) & (dat$R == r),]
  dat_cur <- dat_cur[dat_cur$wtkg == wtkg,]
  dat_cur_Y <- dat_cur[dat_cur$Y == y,]
  res <- dim(dat_cur_Y)[1]/dim(dat_cur)[1]
  if(is.na(res)){
    res <- 0
    print("NaN Pr_Y_SACZR")
  }
  return(res)
}

Pr_R_SACZ <- function(r,s,a,wtkg,z,dat){
  dat_cur <- dat[(dat$S == s) & (dat$A == a) & (dat$Z == z),]
  dat_cur <- dat_cur[dat_cur$wtkg == wtkg,]
  dat_cur_R <- dat_cur[dat_cur$R == r,]
  res <- dim(dat_cur_R)[1]/dim(dat_cur)[1]
  if(is.na(res)){
    res <- 0
    print("NaN: Pr_R_SACZ")
  }
  return(res)
}

int_S_Z <- function(s,z,dat){
  ## compute \int_{A,C} Pr{S = s | Z = z, A, C}f(A,C) d\mu(A)d\mu(C)
  integration <- 0
  for(a in 0:1){
    for(wtkg in 0:1){
        integration = integration + Pr_S_ACZ(s,a,wtkg,z,dat)*f(a,wtkg,dat)
    }
  }
  return(integration)
}

int_Y_SZR <- function(y,s,z,r,dat){
  ## compute \int_{A,C} Pr(Y = y | S = s, Z = z, R =r ,A, C)f(A,C)d\mu(A)d\mu(C)
  integration <- 0
  for(a in 0:1){
    for(wtkg in 0:1){
        integration = integration + Pr_Y_SACZR(y,s,a,wtkg,z,r,dat)*f(a,wtkg,dat)
      }
  }
  return(integration)
}

f_AC_LL <- function(a,wtkg,dat){
  ## compute f(A,C | G = LL)
  res <- Pr_S_ACZ(1,a,wtkg,0,dat)*f(a,wtkg,dat)/int_S_Z(1,0,dat)
  return(res)
}

f_AC_S1 <- function(a,wtkg,dat){
  ## compute f(A,C | S(1) = 1)
  res <- Pr_S_ACZ(1,a,wtkg,1,dat)*f(a,wtkg,dat)/int_S_Z(1,1,dat)
  return(res)
}
##### Compute values
# gamma_proportion <- int_S_Z(1,0,dat)/(1 - int_S_Z(0,1,dat))
gamma_proportion <- function(a,wtkg,dat){
  res <- Pr_S_ACZ(1,a,wtkg,0,dat)/Pr_S_ACZ(1,a,wtkg,1,dat)
  if(is.na(res)){
    res <- 0
    print("gamma proportion: nan")
  }
  if(is.infinite(res)){
    res <- 0
    print("Inf: gamma proportion")
  }
  return(res)
}
for(a in 0:1){
  for(wtkg in 0:1){
      print(gamma_proportion(a,wtkg,dat))
    }
}
for(s in 0:1){
  for(z in 0:1){
    print(int_S_Z(s,z,dat))
  }
}
## a_u
int_a_u <- 0
for(a in 0:1){
  for(wtkg in 0:1){
      pi_x_u <- (Pr_Y_SACZR(1,1,a,wtkg,1,1,dat)*Pr_R_SACZ(1,1,a,wtkg,1,dat)+Pr_R_SACZ(0,1,a,wtkg,1,dat))
      x_gamma <- gamma_proportion(a,wtkg,dat)
      x_a_u <- pi_x_u/x_gamma 
      int_a_u = int_a_u + min(1, x_a_u)*f_AC_LL(a,wtkg,dat)
  }
}


## a_l
int_a_l <- 0
for(a in 0:1){
  for(wtkg in 0:1){
      x_gamma <- gamma_proportion(a,wtkg,dat)
      pi_x_l <- Pr_Y_SACZR(1,1,a,wtkg,1,1,dat)*Pr_R_SACZ(1,1,a,wtkg,1,dat)            
      x_a_l <- pi_x_l/x_gamma - (1 - x_gamma)/x_gamma
      int_a_l = int_a_l + max(0,x_a_l)*f_AC_LL(a,wtkg,dat)
  }
}

int_a_l <- max(0,int_a_l)



#b_u
int_b_u <- 0
for(a in 0:1){
  for(wtkg in 0:1){
      int_b_u <- int_b_u + (Pr_Y_SACZR(1,1,a,wtkg,0,1,dat)*Pr_R_SACZ(1,1,a,wtkg,0,dat) + Pr_R_SACZ(0,1,a,wtkg,0,dat))*f_AC_LL(a,wtkg,dat)
  }
}

int_b_l <- 0
for(a in 0:1){
  for(wtkg in 0:1){
      int_b_l <- int_b_l + (Pr_Y_SACZR(1,1,a,wtkg,0,1,dat)*Pr_R_SACZ(1,1,a,wtkg,0,dat))*f_AC_LL(a,wtkg,dat)
  }
}
print(paste("CE falls in interval: [",int_a_l - int_b_u,int_a_u - int_b_l,"]"))

