# EC Estimand Project -- simulation use dcc
#
# Peijin Wang
# 2025-6-9

library(doParallel)
library(foreach)

rm(list=ls())
setwd("/Users/wangpeijin/Desktop/EC Estimator/Simulation/20250512")

cluster <- 1
nsim <- 1
cl <- makeCluster(cluster)
registerDoParallel(cl)
getDoParWorkers()


## To use argument from BATCH job
args=(commandArgs(TRUE))

# Simulation ####
# parameters
alloc_ratio <- as.numeric(args[1]) # rct allocation ratio
N2 <- as.numeric(args[2])        # ec size
job_id <- as.numeric(args[3])
name <- paste0(args[4], "_", job_id)

alloc_ratio <- 1
N2 <- 100
job_id <- 1


OUT <- foreach(i = 1:cluster, .combine = rbind, .packages = c("tidyverse","PSweight")) %dopar%{
  source("trial_function.R")
  
  start_seed <- (job_id - 1) * (nsim * cluster) + nsim * (i - 1) + 1
  end_seed <- (job_id - 1) * (nsim * cluster) + nsim * (i - 1) + nsim
  
  N1 <- 200
  SIM <- NULL
  for(iter in (nsim*(i-1)+1):(nsim*(i-1)+nsim)){
    for(phi in c(0,0.25,0.5)){
      for(mu.x2.ec in c(0,0.5,1,2)){
        for(sigma.x2.ec in c(sqrt(1),sqrt(1.5))){ 
          N11 <- ifelse(alloc_ratio == 0, N1, N1*alloc_ratio/(alloc_ratio+1))
          SIM <- rbind(SIM, sim_fun(alloc_ratio = alloc_ratio, phi1 = phi, phi2 = phi, mu.x2.ec = mu.x2.ec, sigma.x2.ec = sigma.x2.ec, N2 = N2) %>% 
                         mutate(N1 = N1, N11 = N11, N10 = N1-N11, N2 = N2, lambda = N1/(N1+N2), phi = phi, 
                                mu.x2.ec = mu.x2.ec, sigma.x2.ec = sigma.x2.ec, iter = iter))
        }
      }
    }

  }
  SIM
}


write.csv(OUT, file = paste0("results/OUT_",name,".csv"))



