rm(list = ls())
library(dplyr)
library(gmm)
library(speff2trial)
data("ACTG175")
dat <- ACTG175
Z <- (ACTG175$arms == 1 | ACTG175$arms == 2)
S <- ((ACTG175$cens==0) & (ACTG175$days>7*96))
S[S] <- 1

dat$S <- S
dat_S <- dat[dat$S==1,]

nrow(dat_S)/nrow(dat)
nrow(dat_S[dat_S$r==1,])/nrow(dat_S)

dat_t <- dat[dat$treat==1,]
dat_c <- dat[dat$treat==0,]
dat_tS <- dat_t[dat_t$S==1,]
nrow(dat_tS)/nrow(dat_t)


dat_cS <- dat_c[dat_c$S==1,]
nrow(dat_cS)/nrow(dat_c)

nrow(dat_tS[dat_tS$r==0,])/nrow(dat_tS)
nrow(dat_cS[dat_cS$r==0,])/nrow(dat_cS)

### compute difference
dat_t <- dat_t[dat_t$r == 1,]
dat_c <- dat_c[dat_c$r == 1,]
dat_t <- dat_t[dat_t$S == 1,]
dat_c <- dat_c[dat_c$S == 1,]

dat_t$delta <- ifelse(dat_t$cd496 - dat_t$cd40 > 0,1,0)
dat_c$delta <- ifelse(dat_c$cd496 - dat_c$cd40 > 0,1,0)

sum(dat_t$delta)/nrow(dat_t)
sum(dat_c$delta)/nrow(dat_c)

dat_1 <- dat_t$delta
dat_2 <- dat_c$delta
result <- t.test(dat_1,dat_2)
mean_difference <- result$estimate
md <- mean_difference[1] - mean_difference[2]
conf_interval <- result$conf.int
cat("Mean Difference:", mean_difference, "\n")
cat("95% Confidence Interval for the Difference:", conf_interval[1], "to", conf_interval[2], "\n")

diff <-mean(dat_t$delta) - mean(dat_c$delta)
diff



sum(dat$r)/nrow(ACTG175)

nrow(dat[dat$S==1,])/nrow(dat)

nrow(dat_S[dat_S$r==1,])/nrow(dat_S)



