
# 0. Ideas ----------------------------------------------------------------

# Step 1: process variables(data about 1996)
# Covariates C(pre-treatment variables): age, wtkg, hemo, homo, drugs,
    #karnof, oprior, z30, zprior, preanti, race, gender, str2, strat, symptom
# Note: maybe oprior,z30, zprior can be unioned, and str2, strat

# Shadow variable A: cd40
# Treatment Z: treat
# Outcome Y: cd496

# Survival status: a combination of cens and days.
# If days < 7*96 & cens == 1 , then S = 1. (A little understanding of the definition of S, S == 1 does not means death, but cd4 level too low to satisfy our goal, we want to learn the treatment to normal level cd4, but a low level cd4). Another word, cd4 too low == nearly death(too coarse).
# Else S = 1, but there is an exception, that is when cens==1 & days >= 7*49 = 343(in fact, cd496 is 49 -+ 6 weeks, here, for simplicity, we take 49 weeks)
# So, even there is a cd496 value, the sample may be considered as censoring.

# Missing indicator: if S = 0 and R = 1, missing.

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
# Note 1: covariates could change another choice, especially the union of treatment history.
# case 1: baseline, failed this version: CE = 0.060
# covariates <- c("age","wtkg","hemo","homo","drugs","karnof","gender","str2","symptom")
# case_trial <- 1
# 
# case 2: only age: 0.031
# covariates <- c("age")
# case_trial <- 2
# # case 3: trial: seems okay; ce = 0.014
# covariates <- c("age","wtkg")
# case_trial <- 3
# # case 4: trial: ce = 0.011
# covariates <- c("age","wtkg","karnof","gender","str2","symptom")
# case_trial <- 4
# # case 5: trial: ce = -0.069
# covariates <- c("age","wtkg","karnof","gender")
# case_trial <- 5
# 
# # case 6: trial: ce = -0.030
# covariates <- c("age","wtkg","gender")
# case_trial <- 6
# 
# 
# # case 7: trial: seems okay; ce = 0
# 
# # case 8: trial: ce = 0.081
# covariates <- c("wtkg","gender")
# case_trial <- 8
# 
# # case 9: trial: ce = 0.017
# covariates <- c("wtkg","gender","karnof")
# case_trial <- 9
# 
# # case 10: trial: ce = -0.271
# covariates <- c("wtkg","gender","karnof","str2")
# case_trial <- 10
# 
# case 11: trial: ce = 0.024
covariates <- c("wtkg","gender","karnof","str2","symptom")
case_trial <- 11


# covariates <- c()

C <- rep(1,nrow(ACTG175))
for(c in covariates){
  print(c)
  if(c == 'age' | c == 'wtkg' | c == 'karnof'){
    C <- cbind(C,scale(select(ACTG175,c)))
  }
  else{
    C <- cbind(C,select(ACTG175,c))
  }
}
C <- as.data.frame(C)

colnames(C)[1] <- "bias"
A <- ACTG175$cd40
A <- scale(A)
# Z <- ACTG175$treat


Z <- (ACTG175$arms == 1 | ACTG175$arms == 2)
S <- ((ACTG175$cens==0) & (ACTG175$days>7*96))
S[S] <- 1
# S <- ACTG175$offtrt==0
delta_cd4 <- ACTG175$cd496 - ACTG175$cd40
R <- (ACTG175$r == 1)
R[R] <- 1

medi_delta <- median(delta_cd4,na.rm = TRUE)
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
print("Descriptive Analysis: ")
print(paste("Number of samples:",nrow(ACTG175)))
print(paste("Censoring Proportion:",sum(S)/nrow(ACTG175)))
print(paste("Missing Proportion:",sum(R)/nrow(ACTG175)))

# 3. Model Building -------------------------------------------------------

dat_obs <- dat
### A here, is a problem
## Normalization

dim_status <- c(ncol(C),1,1,1,1,1)



### pr(R = 1 | S = 1 ,A, C, Y; alpha)
### pr(S = 1 | Z = 1, A ,C; beta1)
### pr(S = 1 | Z = 0, A ,C; beta01)
### E[Y | Z = 1, G = g, C ] = \mu_g(C; \gamma_g)
### E[Y | Z = 0, S = 1, C] = m(C; \eta)
## initial values
set.seed(20230701)

ini_val_alpha <- runif( (ncol(C)+2))
ini_val_beta1 <- runif( ncol(C) + 1)
ini_val_beta01 <- runif( ncol(C) + 1)
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
    data = select(dat_obs, c('bias',covariates,'A','Z','S')),
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

### singular:: to check the gradients. rank(A*A^T) = rank(A), A is gradients of m_sum
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

# 4. Estimation of causal effect ------------------------------------------
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
CE <- E1LL - E0LL
print(CE)


# 
# 
# # 5. Boostrap Compute Variance --------------------------------------------
library(parallel)
clnum <- detectCores() - 1
print(paste("detectCores:" ,clnum))
print(paste("specified cores:" , clnum))
print(clnum)



recycle_times <- 200
seeds <- 1
# 
# 
# 
# 
# 
# 5. Bootstrap to compute quantiles ---------------------------------------
cl <- makeCluster(getOption('cl.cores',clnum))
print("makecluster starts")
time_start <- Sys.time()


res <- mcmapply(Estimation_boot,seeds = 1:recycle_times,MoreArgs = list(dat_obs = dat_obs,dim_status = dim_status),mc.cores = clnum)
time_end <- Sys.time()
print(time_end - time_start)
CEs_bootstrap <- c()
for(i in 1:length(res)){
  CEs_bootstrap <- c(CEs_bootstrap,res[i])
}
stopCluster(cl)
save.image(paste0("Estimation_",case_trial,".RData"))
print(paste0(case_trial, " is over."))


