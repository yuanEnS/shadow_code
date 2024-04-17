library(speff2trial)
data("ACTG175")
head(ACTG175)
help("ACTG175")
# 
# A data frame with 2139 observations on the following 27 variables:
#   
# pidnum: patient's ID number
# age: age in years at baseline
# wtkg: weight in kg at baseline
# hemo: hemophilia (0=no, 1=yes)
# homo: homosexual activity (0=no, 1=yes)
# drugs: history of intravenous drug use (0=no, 1=yes)
# karnof: Karnofsky score (on a scale of 0-100)
# oprior: non-zidovudine antiretroviral therapy prior to initiation of study treatment (0=no, 1=yes)
# z30: zidovudine use in the 30 days prior to treatment initiation (0=no, 1=yes)
# zprior: zidovudine use prior to treatment initiation (0=no, 1=yes)
# preanti: number of days of previously received antiretroviral therapy
# race: race (0=white, 1=non-white)
# gender: gender (0=female, 1=male)
# str2: antiretroviral history (0=naive, 1=experienced)
# strat: antiretroviral history stratification (1='antiretroviral naive', 2='> 1 but ≤≤ 52 weeks of prior antiretroviral therapy', 3='> 52 weeks')
# symptom: symptomatic indicator (0=asymptomatic, 1=symptomatic)

# treat: treatment indicator (0=zidovudine only, 1=other therapies)
# offtrt: indicator of off-treatment before 96±±5 weeks (0=no,1=yes)
# cd40: CD4 T cell count at baseline
# cd420: CD4 T cell count at 20±±5 weeks
# cd496: CD4 T cell count at 96±±5 weeks (=NA if missing)
# r: missing CD4 T cell count at 96±±5 weeks (0=missing, 1=observed)
# cd80: CD8 T cell count at baseline
# cd820: CD8 T cell count at 20±±5 weeks
# cens: indicator of observing the event in days
# days: number of days until the first occurrence of: (i) a decline in CD4 T cell count of at least 50 (ii) an event indicating progression to AIDS, or (iii) death.
# arms: treatment arm (0=zidovudine, 1=zidovudine and didanosine, 2=zidovudine and zalcitabine, 3=didanosine).

dat <- ACTG175
sum(is.na(dat$cd40))
head(dat)

## Why cd496
attach(ACTG175)
### reason 1: time should not be too long due to the limit of study resources, 2 years, the cd4 level will get stable, generally. 
# reference: https://i-base.info/guides/art-in-pictures/hiv-after-starting-art


### reason 2: if possible, more data, more better, if time is too long, there will be more missing data, which will affect the efficiency of our study conclusions, and censoring is also a problem
### reason 3: longer time, more resources.
hist(days) ## after 96 +- 6 weeks, there are many censoring, to save information, the time is better before then.
dat_treat <- dat[dat$treat == 1,]
dat_control <- dat[dat$treat == 0,]
hist(dat_treat$days)
hist(dat_control$days)
hist(cd496)
hist(dat_treat$cd496)
hist(dat_control$cd496)
