
# 1. N = 2000 -------------------------------------------------------------
# N = 2000, cyc times = 200 Done!
load('2_Simulations/Repeat200_N2000.RData')
true_causal <- 0.3322807
bias <- mean(CEs) - true_causal
bias
rmse <- sqrt(mean((CEs - true_causal)^2))
rmse
coverage_rate <- sum(coverages)/200
coverage_rate

# N = 2000, cyc times = 500 Done!
load('2_Simulations/Repeat500_N2000.RData')
cyc <- 500
true_causal <- 0.3322807
CEs <- NULL
E1ls <- NULL
E2ls <- NULL
gamma_LL1 <- NULL
hat_eta_1s <- NULL
hat_eta_2s <- NULL
coverages <- NULL

for(i in 1:(cyc+error_trial)){
  if(length(CEs) == cyc){
    break
  }
  tryCatch(
    expr = {
      CEs <- c(CEs,res[[i]][1]-res[[i]][2])
      E1ls <- c(E1ls,res[[i]][1])
      E2ls <- c(E2ls,res[[i]][2])
      gamma_LL1 <- c(gamma_LL1,res[[i]][3])
      hat_eta_1s <- c(hat_eta_1s,res[[i]][7])
      hat_eta_2s <- c(hat_eta_2s,res[[i]][8])
      coverages <- c(coverages,res[[i]][19])
      message(paste("Successfully added cyc: ",i))
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
    },
    finally = {
      message('All done, quitting.')
    }
  )
}

bias <- mean(CEs) - true_causal
rmse <- sqrt(mean((CEs - true_causal)^2))
coverage_rate <- sum(coverages)/cyc
bias
rmse
coverage_rate

# N = 2000, cyc times = 1000 Done!
load('2_Simulations/Repeat1000_N2000.RData')
cyc <- 1000
true_causal <- 0.3322807
CEs <- NULL
E1ls <- NULL
E2ls <- NULL
gamma_LL1 <- NULL
hat_eta_1s <- NULL
hat_eta_2s <- NULL
coverages <- NULL

for(i in 1:(cyc+error_trial)){
  if(length(CEs) == cyc){
    break
  }
  tryCatch(
    expr = {
      CEs <- c(CEs,res[[i]][1]-res[[i]][2])
      E1ls <- c(E1ls,res[[i]][1])
      E2ls <- c(E2ls,res[[i]][2])
      gamma_LL1 <- c(gamma_LL1,res[[i]][3])
      hat_eta_1s <- c(hat_eta_1s,res[[i]][7])
      hat_eta_2s <- c(hat_eta_2s,res[[i]][8])
      coverages <- c(coverages,res[[i]][19])
      message(paste("Successfully added cyc: ",i))
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
    },
    finally = {
      message('All done, quitting.')
    }
  )
}

bias <- mean(CEs) - true_causal
rmse <- sqrt(mean((CEs - true_causal)^2))
coverage_rate <- sum(coverages)/cyc
bias
rmse
coverage_rate


# 2. N = 5000 -------------------------------------------------------------
# collect data
rm(list = ls())
CE_5000 <- c()
coverages_5000 <- c()
load('2_Simulations/Repeat200_N5000.RData')
CE_5000 <- c(CE_5000,CEs) 
coverages_5000 <- c(coverages_5000,coverages)
load('2_Simulations/Repeat200_N5000_1.RData')
CE_5000 <- c(CE_5000,CEs) 
coverages_5000 <- c(coverages_5000,coverages)
load('2_Simulations/Repeat200_N5000_2.RData')
CE_5000 <- c(CE_5000,CEs) 
coverages_5000 <- c(coverages_5000,coverages)
load('2_Simulations/Repeat200_N5000_3.RData')
CE_5000 <- c(CE_5000,CEs) 
coverages_5000 <- c(coverages_5000,coverages)
load('2_Simulations/Repeat200_N5000_4.RData')
CE_5000 <- c(CE_5000,CEs) 
coverages_5000 <- c(coverages_5000,coverages)


# cyc times = 200
cyc_times <- 200
CEs_200 <- CE_5000[1:cyc_times]
coverages_200 <- coverages_5000[1:cyc_times]
true_causal <- 0.3322807
CEs <- CEs_200
coverages <- coverages_200

bias <- mean(CEs) - true_causal
bias
rmse <- sqrt(mean((CEs - true_causal)^2))
rmse
coverage_rate <- sum(coverages)/200
coverage_rate

# cyc times = 500
cyc_times <- 500
CEs <- CE_5000[1:cyc_times]
coverages <- coverages_5000[1:cyc_times]
true_causal <- 0.3322807

bias <- mean(CEs) - true_causal
bias
rmse <- sqrt(mean((CEs - true_causal)^2))
rmse
coverage_rate <- sum(coverages)/cyc_times
coverage_rate

# cyc times = 1000
cyc_times <- 1000
CEs <- CE_5000[1:cyc_times]
coverages <- coverages_5000[1:cyc_times]
true_causal <- 0.3322807

bias <- mean(CEs) - true_causal
bias
rmse <- sqrt(mean((CEs - true_causal)^2))
rmse
coverage_rate <- sum(coverages)/cyc_times
coverage_rate

# 3. N = 8000 -------------------------------------------------------------
# collect data
rm(list = ls())
CE_8000 <- c()
coverages_8000 <- c()
load('2_Simulations/Repeat200_N8000.RData')
CE_8000 <- c(CE_8000,CEs) 
coverages_8000 <- c(coverages_8000,coverages)
load('2_Simulations/Repeat200_N8000_1.RData')
CE_8000 <- c(CE_8000,CEs) 
coverages_8000 <- c(coverages_8000,coverages)
load('2_Simulations/Repeat200_N8000_2.RData')
CE_8000 <- c(CE_8000,CEs) 
coverages_8000 <- c(coverages_8000,coverages)
load('2_Simulations/Repeat200_N8000_3.RData')
CE_8000 <- c(CE_8000,CEs) 
coverages_8000 <- c(coverages_8000,coverages)
load('2_Simulations/Repeat200_N8000_4.RData')
CE_8000 <- c(CE_8000,CEs) 
coverages_8000 <- c(coverages_8000,coverages)


# cyc times = 200
cyc_times <- 200
CEs_200 <- CE_8000[1:cyc_times]
coverages_200 <- coverages_8000[1:cyc_times]
true_causal <- 0.3322807
CEs <- CEs_200
coverages <- coverages_200

bias <- mean(CEs) - true_causal
bias
rmse <- sqrt(mean((CEs - true_causal)^2))
rmse
coverage_rate <- sum(coverages)/cyc_times
coverage_rate

# cyc times = 500
cyc_times <- 500
CEs <- CE_8000[1:cyc_times]
coverages <- coverages_8000[1:cyc_times]
true_causal <- 0.3322807

bias <- mean(CEs) - true_causal
bias
rmse <- sqrt(mean((CEs - true_causal)^2))
rmse
coverage_rate <- sum(coverages)/cyc_times
coverage_rate


# cyc times = 1000
cyc_times <- 1000
CEs <- CE_8000[1:cyc_times]
coverages <- coverages_8000[1:cyc_times]
true_causal <- 0.3322807

bias <- mean(CEs) - true_causal
bias
rmse <- sqrt(mean((CEs - true_causal)^2))
rmse
coverage_rate <- sum(coverages)/cyc_times
coverage_rate




# 4. Bound Results --------------------------------------------------------

### Bound
load("2_Simulations/Bound/Repeat_bound_ori_gene1000_N8000.RData")
lower_bounds <- NULL
upper_bounds <- NULL
len_intervals <- NULL
CEs <- NULL
for(i in 1:(cyc+error_trial)){
  if(length(lower_bounds) == 1000){
    break
  }
  tryCatch(
    expr = {
      lower_bounds <- c(lower_bounds,res[[i]][1])
      upper_bounds <- c(upper_bounds,res[[i]][2])
      len_intervals <- c(len_intervals,res[[i]][3])
      CEs <- c(CEs,res[[i]][4])
      message(paste("Successfully added cyc: ",i))
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
    },
    finally = {
      message('All done, quitting.')
    }
  )
}
length(CEs)
hist(lower_bounds)
hist(upper_bounds)
hist(len_intervals)

true_causal <- 0.3322807
data <- data.frame(Experiment_number=rep(1:length(lower_bounds), 3),
                   Causal_Effect=c(upper_bounds,CEs,lower_bounds), 
                   variable=c(rep("Upper bound",length(lower_bounds)),rep("Estimated causal effect",length(lower_bounds)),rep("Lower bound",length(lower_bounds))))




library(ggplot2)
# Create the violin plot

ggplot(data, aes(x = variable, y = Causal_Effect)) +
  geom_violin(aes(fill = variable)) +  # Violin plot
  geom_hline(yintercept = true_causal, color = "black", linetype = "dashed") +  # Add horizontal line
  labs(x = "Group", y = "Value", title = "Causal Effects") +
  theme_classic()
ggsave("./2_Simulations/Bound/Bounds_violin.png", width = 8)
