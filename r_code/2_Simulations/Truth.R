### compute true E(Y_0 | G = LL ) and E(Y_1 | G = LL )
## simulation setting

library(parallel)

clnum <- detectCores()
cl <- makeCluster(getOption('cl.cores', clnum))

source("./1_functions/M1_data_generating.R")
cyc <- 100
res1 <- parLapply(cl, 1:cyc, data_generate)
dat <- NULL
for (i in 1:cyc) {
  dat <- rbind(dat, res1[[i]])
}
df_LL <- dat[dat$G == 0 ,]
E1_true <- mean(df_LL$Y_1_ll)
E0_true <- mean(df_LL$Y_0_ll)
print(E1_true)
print(E0_true)
print(E1_true -  E0_true)

## truth = 0.3322807