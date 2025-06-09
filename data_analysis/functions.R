# EC Estimand Paper -- case study
# functions
# 
# Peijin Wang
# 2025-6-9


#### Table 1 ####
# add stabilized weights to data
data_sw_fun <- function(trial, ec, type = "PATE"){
  # type takes ("PATE","PATT","PATC","PATO") = (ATI, ATT, ATEC, ATO)
  
  data <- rbind(trial, ec)
  
  # estimate trial propensity score Pr(Z=1|X)
  data$ps <- SumStat(ps.formula = Z ~ AGEREF + DIAGDYS + UVSTRAT + BLINJECT +
                       BLNJNT_1 + BLNJNT_2 + BLNJNT_3 + JADAS10,
                     trtgrp = "1", zname = "Z", data = data)$propensity[,2]
  
  # weights
  if(type == "PATE"){ data$w1 <- 1/data$ps; data$w0 <- 1/(1-data$ps) }
  if(type == "PATT"){ data$w1 <- rep(1,nrow(data)); data$w0 <- data$ps/(1-data$ps) }
  if(type == "PATC"){ data$w1 <- (1-data$ps)/data$ps; data$w0 <- rep(1,nrow(data)) }
  if(type == "PATO"){ data$w1 <- 1-data$ps; data$w0 <- data$ps }
  # stabilized weight
  data$sw1 <- mean(data$Z==1) * data$w1; data$sw0 <- mean(data$Z==0) * data$w0
  
  return(data)
}

# Table 1 of filtered data
## compute original/weighted sample's summary statistics
sw_mean_fun <- function(data, covariate, type = "raw", pre_type = "raw"){
  # type takes ("PATE","PATT","PATC","PATO") = (ATI, ATT, ATEC, ATO)
  # pre_type takes ("raw", "combine", "by Z")
  
  if(pre_type == "raw"){
    results <- cbind(data %>% filter(Z == 1) %>% 
                       summarise(trial_mean = mean(!!sym(covariate)),
                                 trial_q50 = quantile(!!sym(covariate), probs = 0.5),
                                 trial_q25 = quantile(!!sym(covariate), probs = 0.25),
                                 trial_q75 = quantile(!!sym(covariate), probs = 0.75)),
                     data %>% filter(Z == 0) %>% 
                       summarise(EC_mean = mean(!!sym(covariate)),
                                 EC_q50 = quantile(!!sym(covariate), probs = 0.5),
                                 EC_q25 = quantile(!!sym(covariate), probs = 0.25),
                                 EC_q75 = quantile(!!sym(covariate), probs = 0.75))) %>%
      mutate(covariate = covariate)
  }
  if(pre_type == "combine"){ # weighted + stabilized weight (excluded from paper)
    # stabilized weight
    data$norm_sw <- with(data, Z*sw1/sum(Z*w1) + (1-Z)*sw0/sum((1-Z)*w0))
    
    # weighted mean
    results <- data %>% summarise(
      sw_X_mean = sum(norm_sw*!!sym(covariate)),
      sw_X_q50 = wtd.quantile(!!sym(covariate), weights = norm_sw, probs = 0.50, normwt = T),
      sw_X_q25 = wtd.quantile(!!sym(covariate), weights = norm_sw, probs = 0.25, normwt = T),
      sw_X_q75 = wtd.quantile(!!sym(covariate), weights = norm_sw, probs = 0.75, normwt = T)) %>% 
      mutate(covariate = covariate)
    
    colnames(results)[colnames(results) == "sw_X_mean"] <- paste0("sw_X_mean_",type)
    colnames(results)[colnames(results) == "sw_X_q50"] <- paste0("sw_X_q50_",type)
    colnames(results)[colnames(results) == "sw_X_q25"] <- paste0("sw_X_q25_",type)
    colnames(results)[colnames(results) == "sw_X_q75"] <- paste0("sw_X_q75_",type)
    
  }
  if(pre_type == "by Z") { # weighted + standardized weight
    # standardized weight
    data$norm_w <- with(data, Z*w1/sum(Z*w1) + (1-Z)*w0/sum((1-Z)*w0))
    
    # weighted mean
    results <- data %>% group_by(Z) %>% summarise(
      sw_X_mean = sum(norm_w*!!sym(covariate)),
      sw_X_q50 = wtd.quantile(!!sym(covariate), weights = norm_w, probs = 0.50, normwt = T),
      sw_X_q25 = wtd.quantile(!!sym(covariate), weights = norm_w, probs = 0.25, normwt = T),
      sw_X_q75 = wtd.quantile(!!sym(covariate), weights = norm_w, probs = 0.75, normwt = T)) %>% 
      mutate(covariate = covariate) %>% 
      pivot_wider(names_from = Z, values_from = c(sw_X_mean, sw_X_q50, sw_X_q25, sw_X_q75))
    
    colnames(results)[colnames(results) == "sw_X_mean_0"] <- paste0("sw_X_mean_0_",type)
    colnames(results)[colnames(results) == "sw_X_mean_1"] <- paste0("sw_X_mean_1_",type)
    
    colnames(results)[colnames(results) == "sw_X_q50_0"] <- paste0("sw_X_q50_0_",type)
    colnames(results)[colnames(results) == "sw_X_q50_1"] <- paste0("sw_X_q50_1_",type)
    
    colnames(results)[colnames(results) == "sw_X_q25_0"] <- paste0("sw_X_q25_0_",type)
    colnames(results)[colnames(results) == "sw_X_q25_1"] <- paste0("sw_X_q25_1_",type)
    
    colnames(results)[colnames(results) == "sw_X_q75_0"] <- paste0("sw_X_q75_0_",type)
    colnames(results)[colnames(results) == "sw_X_q75_1"] <- paste0("sw_X_q75_1_",type)
  }
  return(results)
}

## generate table 1
table1_fun <- function(rct, ec, covariates, type = "raw", pre_type = "raw"){
  # type takes ("PATE","PATT","PATC","PATO") = (ATI, ATT, ATEC, ATO)
  # pre_type takes ("raw", "combine", "by Z")
  
  if(pre_type == "raw"){ data <- rbind(rct, ec) }
  if(type != "raw"){ data <- data_sw_fun(rct, ec, type) }
  
  map_df(covariates, ~sw_mean_fun(data, .x, type, pre_type))
}





#### EC Estimators ####

# compute estimator
ec_est_fun <- function(rct, ec, type = "PATE"){
  # type takes ("PATE","PATT","PATC","PATO") = (ATI, ATT, ATEC, ATO)

  N10 <- nrow(rct[rct$A==0,])
  N2 <- nrow(ec)
  data <- rbind(rct,ec)
  
  # estimate trial propensity score Pr(Z=1|X)
  data$ps <- SumStat(ps.formula = Z ~ AGEREF + DIAGDYS + UVSTRAT + BLINJECT + 
                       BLNJNT_1 + BLNJNT_2 + BLNJNT_3 + rcs(JADAS10, df),
                     trtgrp = "1", zname = "Z", data = data)$propensity[,2]
  
 
  # weights
  if(type == "PATE"){ data$w1 <- 1/data$ps; data$w0 <- 1/(1-data$ps) }
  if(type == "PATT"){ data$w1 <- rep(1,nrow(data)); data$w0 <- data$ps/(1-data$ps) }
  if(type == "PATC"){ data$w1 <- (1-data$ps)/data$ps; data$w0 <- rep(1,nrow(data)) }
  if(type == "PATO"){ data$w1 <- 1-data$ps; data$w0 <- data$ps }
  
  # estimator
  estimator <- with(data, sum(w1*A*Z*Y)/sum(w1*A*Z) - (N10/(N10+N2))*sum(w1*(1-A)*Z*Y)/sum(w1*(1-A)*Z) 
                    - (N2/(N10+N2))*sum(w0*(1-A)*(1-Z)*Y)/sum(w0*(1-A)*(1-Z))  )
  # normalized weights
  data <- data %>% mutate(nw = case_when(Z==1 & A==1 ~ w1/sum(w1[Z==1&A==1]),
                                         Z==1 & A==0 ~ (N10/(N10+N2))*w1/sum(w1[Z==1&A==0]),
                                         Z==0        ~ (N2/(N10+N2))*w0/sum(w0[Z==0]),
                                         TRUE ~ NA_real_)) %>% mutate(Y.weighted = nw*Y)

  return(list(data = data,estimator = estimator))
}

# compute standard error and CI
ec_ci_fun <- function(rct, ec, type = "PATE", n.boot = 1000){
  # bootstrap computing se
  ## point estimate
  point.est <- ec_est_fun(rct,ec,type)[["estimator"]]

  n.rct <- nrow(rct)
  n.ec <- nrow(ec)
  boot.est <- replicate(n.boot,{
    rct.id.boot <- sample(1:n.rct, n.rct, replace = T)
    ec.id.boot <- sample(1:n.ec, n.ec, replace = T)
    ec_est_fun(rct[rct.id.boot,],ec[ec.id.boot,],type)[["estimator"]]
  })
  boot.se <- sd(boot.est)
  boot.lb <- quantile(boot.est, probs = 0.025)
  boot.ub <- quantile(boot.est, probs = 0.975)
  return(tibble(est = point.est, se = boot.se, lower = boot.lb, upper = boot.ub))
}


#### Back up functions (not used in the paper) #######

#### RCT Estimator ####
rct_est_fun <- function(rct, model){
  # model = c("interaction","no interaction","non-parametric")
  # center 
  rct_centered <- rct %>%
    mutate(across(where(is.numeric) & !c(Y, Z, A), scale, scale = FALSE))
  
  if(model == "non-parametric"){ fit <- lm(Y ~ A, data = rct_centered) }
  if(model == "no interaction"){
    fit <- lm(Y ~ A + age + female + whit + duripf + bfvcpp + bsdlcopp + reflux + cadhist + bsfev1fvc, 
              data = rct_centered)
  }
  if(model == "interaction"){
    fit <- lm(Y ~ A * (age + female + whit + duripf + bfvcpp + bsdlcopp + reflux + cadhist + bsfev1fvc), 
              data = rct_centered)
  }
  # Extract the coefficient of interest
  estimator <- coef(fit)[2]
  return(estimator)
}

# compute standard error and CI
rct_ci_fun <- function(rct, interaction = T, n.boot = 1000){
  # bootstrap computing se
  ## point estimate
  point.est <- rct_est_fun(rct, interaction)
  
  n.rct <- nrow(rct)
  boot.est <- replicate(n.boot,{
    rct.id.boot <- sample(1:n.rct, n.rct, replace = T)
    rct_est_fun(rct[rct.id.boot,], interaction)
  })
  boot.se <- sd(boot.est)
  boot.lb <- quantile(boot.est, probs = 0.025)
  boot.ub <- quantile(boot.est, probs = 0.975)
  return(tibble(est = point.est, se = boot.se, lower = boot.lb, upper = boot.ub))
}