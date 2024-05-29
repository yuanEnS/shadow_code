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


source("./r_code//functions4data/MLEs.R")

covariates <- c("wtkg",'karnof')


median_wtkg <- median(select(ACTG175,'wtkg')$wtkg)
C <- ifelse((select(ACTG175,'wtkg'))$wtkg>median_wtkg,1,0)
C <- cbind(C,ifelse((select(ACTG175,'karnof'))$karnof>95,1,0))

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
C <- as.data.frame(C)
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

delta_z <- function(z,dat){
  dat_cur <- dat[dat$S == 1,]
  dat_cur <- dat_cur[dat_cur$Z == z,]
  dat_final <- dat_cur[dat_cur$R ==1, ]
  res <- dim(dat_final)[1]/dim(dat_cur)[1]
  return(res)
}
delta_1 <- delta_z(1,dat)

delta_0 <- delta_z(0,dat)

gamma_com <- function(dat){
  dat_z_0 <- dat[dat$Z==0,]
  
  dat_z_1 <- dat[dat$Z==1,]
  
  dat_cur_0 <- dat_z_0[dat_z_0$S==1,]
  dat_cur_1 <- dat_z_1[dat_z_1$S==1,]
  res <- (dim(dat_cur_0)[1]/dim(dat_z_0)[1])/((dim(dat_cur_1)[1]/dim(dat_z_1)[1]))
  return(res)
}

gamma <- gamma_com(dat)

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

print(paste("CE falls in interval: [",theta_111_l - theta_011_u,theta_111_u - theta_011_l,"]"))

