
# 0. Ideas ----------------------------------------------------------------

# Maybe the true causal effect is just 0.

# 1. Packages -------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(gmm)
library(speff2trial)

# 2. Data processing ------------------------------------------------------
data("ACTG175")
attach(ACTG175)

# Corase analysis
print("coarse analysis:")
print(paste("Mean cd4 delta of treatment group is:",mean(ACTG175[treat==1,]$cd496 - ACTG175[treat==1,]$cd40,na.rm = TRUE)))
print(paste("Mean cd4 delta of control group is:",mean(ACTG175[treat==0,]$cd496 - ACTG175[treat==0,]$cd40,na.rm = TRUE)))

print("Consider censoring, that is, missing is MCAR: ")
S <- (ACTG175$cens & (ACTG175$days<7*49))
S[S] <- 1
print(paste("Mean cd496 of treatment group(consider censoring) is:",mean(ACTG175[treat==1 & S == 1,]$cd496 - ACTG175[treat==1 & S == 1,]$cd40,na.rm = TRUE)))
print(paste("Mean cd496 of control group(consider censoring) is:",mean(ACTG175[treat==0 & S == 1,]$cd496 - ACTG175[treat==0 & S == 1,]$cd40,na.rm = TRUE)))

sum(S)/nrow(ACTG175)

