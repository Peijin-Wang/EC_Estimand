# EC Estimand Project -- functions
#
# Peijin Wang
# 2025-6-9


#### Generate Super Population -- for estimands ####
# generate trial and EC super population for estimand computation
# return super population
super_pop_fun <- function(N1 = 10^6, lambda = 0.5, b0 = 0, beta = c(1, 1), 
                          phi1 = 0, phi2 = 0, p.x1 = 0.5, mu.x2 = 0, sigma.x2 = 1, b.trt = 0, sigma = 1, 
                          source = "trial") {
  # N1: trial super population size; lambda: mixture proportion used to compute ec population size
  # source: ("trial","ec") if trial no need for lambda
  # b0, beta, phi1, phi2: potential outcome model coefficient
  # p.x1, mu.x2, sigma.x2, sigma: covariates and random error dist parameters
  # b.trt: treatment effect

  # EC super population size
  if (source == "ec") { N <- (1 - lambda) / lambda * N1 }
  if (source == "trial") { N <- N1 }

  X1 <- rbinom(N, 1, p.x1)
  X2 <- rnorm(N, mu.x2, sigma.x2)
  b1 <- beta[1]
  b2 <- beta[2]

  # potential outcomes
  Y0 <- b0 + b1 * X1 + b2 * X2 + rnorm(N, 0, sigma)
  Y1 <- Y0 + b.trt + phi1 * X1 + phi2 * X2

  data <- data.frame(Y0 = Y0, Y1 = Y1, X1 = X1, X2 = X2, source = source)
  return(data)
}

# ~ Add true propensity score ####
# compute true propensity score for oberved data
# return full observed data with true PS
pi_fun <- function(pop_trial, pop_ec, p.x1 = c(0.5, 0.5), mu.x2 = c(0, 0), sigma.x2 = c(1, 1), lambda = 0.5) {
  # pop_trial: observed trial
  # pop_ec: observed ec
  # p.x1: X1 parameter for trial and ec
  # mu.x2: X2 mean parameter for trial and ec
  # sigma.x2: X2 sd parameter for trial and ec
  # lambda: mixture proportion

  data <- rbind(pop_trial, pop_ec) %>%
    mutate(
      f1.x1 = dbinom(X1, 1, p.x1[1]), f1.x2 = dnorm(X2, mu.x2[1], sigma.x2[1]),
      f2.x1 = dbinom(X1, 1, p.x1[2]), f2.x2 = dnorm(X2, mu.x2[2], sigma.x2[2]),
      f1 = f1.x1 * f1.x2, f2 = f2.x1 * f2.x2,
      pi = lambda * f1 / (lambda * f1 + (1 - lambda) * f2)
    )
  return(data)
}


#### Estimand ####
# compute estimands
# return estimand and super population used to compute estimand
estimand_fun <- function(super_pop, type = "PATE") {
  # super population (with pi)
  # type = c("PATE","PATT","PATEC", "PATO") corresponding to ATI, ATT, ATEC and ATO in manuscript

  # tilting function
  super_pop$h <- case_when(type == "PATE"  ~ 1,
                           type == "PATT"  ~ super_pop$pi,
                           type == "PATEC" ~ 1 - super_pop$pi,
                           type == "PATO"  ~ super_pop$pi * (1 - super_pop$pi))
  # estimand
  super_pop$cate <- super_pop$Y1 - super_pop$Y0
  estimand <- sum(super_pop$h * super_pop$cate) / sum(super_pop$h)
  
  return(list(estimand = estimand, super_pop = super_pop))
}


#### Generate samples ####
sample_fun <- function(N, alloc_ratio = 1, b0 = 0, beta = c(1, 1), phi1 = 0, phi2 = 0, p.x1 = 0.5, mu.x2 = 0, sigma.x2 = 1, b.trt = 0, sigma = 1, source = "trial") {
  # N: total sample size of trial (N1 in manuscript), of ec (N2 in manuscript)
  # alloc_ratio: trial allocation ratio
  # source = c("trial","ec")

  # generate observed samples
  X1 <- rbinom(N, 1, p.x1)
  X2 <- rnorm(N, mu.x2, sigma.x2)
  b1 <- beta[1]
  b2 <- beta[2]

  # potential outcomes
  Y0 <- b0 + b1 * X1 + b2 * X2 + rnorm(N, 0, sigma)
  Y1 <- Y0 + b.trt + phi1 * X1 + phi2 * X2

  data <- data.frame(Y0 = Y0, Y1 = Y1, X1 = X1, X2 = X2, source = source) %>% mutate(Z = ifelse(source == "trial", 1, 0))

  # observed Y
  if (source == "trial") {
    if(alloc_ratio == 0) { # single-arm
      data$A <- 1
    }else { # rct
      select_A <- sample(1:nrow(data), size = N * alloc_ratio / (alloc_ratio + 1), replace = F)
      data$A <- 0
      data$A[select_A] <- 1
    }
  }

  if (source == "ec") {
    data$A <- 0
  }

  data$Y <- data$A * data$Y1 + (1 - data$A) * data$Y0
  return(data)
}


#### Estimator ####
# compute estimator
ec_est_fun <- function(trial, ec, type = "PATE") {
  # input trial and ec observed data
  # type takes ("PATE","PATT","PATEC","PATO") corresponding to ATI, ATT, ATEC and ATO in manuscript
  
  N10 <- nrow(trial[trial$A == 0, ])
  N2 <- nrow(ec)
  data <- rbind(trial, ec)
  # estimate trial propensity score Pr(Z=1|X)
  data$ps <- SumStat(ps.formula = Z ~ X1 + X2, trtgrp = "1", zname = "Z", data = data)$propensity[, 2]
  # weights
  if (type == "PATE") {
    data$w1 <- 1 / data$ps
    data$w0 <- 1 / (1 - data$ps)
  }
  if (type == "PATT") {
    data$w1 <- rep(1, nrow(data))
    data$w0 <- data$ps / (1 - data$ps)
  }
  if (type == "PATEC") {
    data$w1 <- (1 - data$ps) / data$ps
    data$w0 <- rep(1, nrow(data))
  }
  if (type == "PATO") {
    data$w1 <- 1 - data$ps
    data$w0 <- data$ps
  }
  # estimator
  if(N10 == 0) { # single-arm
    estimator <- with(data, sum(w1 * A * Z * Y) / sum(w1 * A * Z) - sum(w0 * (1 - A) * (1 - Z) * Y) / sum(w0 * (1 - A) * (1 - Z)))
  }else { # rct
    estimator <- with(data, sum(w1 * A * Z * Y) / sum(w1 * A * Z) - (N10 / (N10 + N2)) * sum(w1 * (1 - A) * Z * Y) / sum(w1 * (1 - A) * Z)
                      - (N2 / (N10 + N2)) * sum(w0 * (1 - A) * (1 - Z) * Y) / sum(w0 * (1 - A) * (1 - Z)))
  }
  return(estimator)
}

# compute standard error and CI
# (not included in the simulation)
ec_ci_fun <- function(trial, ec, type = "PATE", n.boot = 1) {
  # bootstrap computing se
  ## point estimate
  point.est <- ec_est_fun(trial, ec, type)

  n.trial <- nrow(trial)
  n.ec <- nrow(ec)
  boot.est <- replicate(n.boot, {
    trial.id.boot <- sample(1:n.trial, n.trial, replace = T)
    ec.id.boot <- sample(1:n.ec, n.ec, replace = T)
    ec_est_fun(trial[trial.id.boot, ], ec[ec.id.boot, ], type)
  })
  boot.se <- sd(boot.est)
  boot.lb <- quantile(boot.est, probs = 0.025)
  boot.ub <- quantile(boot.est, probs = 0.975)
  return(tibble(est = point.est, se = boot.se, lower = boot.lb, upper = boot.ub))
}


#### trial Estimator ####
rct_est_fun <- function(trial) {
  # trial data
  # Extract the coefficient of interest
  fit <- lm(Y ~ A, data = trial)
  estimator <- coef(fit)[2]
  return(estimator)
}

# compute standard error and CI
rct_ci_fun <- function(trial, interaction = T, n.boot = 1) {
  # bootstrap computing se
  ## point estimate
  point.est <- rct_est_fun(trial, interaction)

  n.trial <- nrow(trial)
  boot.est <- replicate(n.boot, {
    trial.id.boot <- sample(1:n.trial, n.trial, replace = T)
    rct_est_fun(trial[trial.id.boot, ], interaction)
  })
  boot.se <- sd(boot.est)
  boot.lb <- quantile(boot.est, probs = 0.025)
  boot.ub <- quantile(boot.est, probs = 0.975)
  return(tibble(est = point.est, se = boot.se, lower = boot.lb, upper = boot.ub))
}



#### Simulation ####
# compute estimand and estimator
output_fun <- function(super_pop, N1 = 200, alloc_ratio = 1, N2 = 100, phi1 = 0, phi2 = 0, mu.x2.ec = 0, sigma.x2.ec = 1, type = "PATE") {
  # super_pop: super population with i
  # N1: sample size of trial, alloc_ratio: trial allocation ratio
  # N2: sample size of ec
  # type: estimand type c("PATE","PATT","PATO)

  # estimand
  estimand <- estimand_fun(super_pop, type)$estimand
  # estimator
  trial <- sample_fun(N = N1, alloc_ratio = alloc_ratio, phi1 = phi1, phi2 = phi2, source = "trial")
  ec <- sample_fun(N = N2, phi1 = phi1, phi2 = phi2, mu.x2 = mu.x2.ec, sigma.x2 = sigma.x2.ec, source = "ec")
  est <- ec_ci_fun(trial, ec, type)

  return(data.frame(type = type, estimand = estimand, estimator = est$est, se = est$se, lower = est$lower, upper = est$upper))
}

# run simulation once
sim_fun <- function(phi1 = 0, phi2 = 0, mu.x2.ec = 0, sigma.x2.ec = 1, b.trt = 0,
                    N1 = 200, alloc_ratio = 1, N2 = 50, pop_only = F) {
  # mixture proportion
  lambda <- N1 / (N1 + N2)
  # generate super population
  pop_rct <- super_pop_fun(phi1 = phi1, phi2 = phi2, b.trt = b.trt, source = "trial")
  pop_ec <- super_pop_fun(lambda = lambda, phi1 = phi1, phi2 = phi2, mu.x2 = mu.x2.ec, sigma.x2 = sigma.x2.ec, b.trt = b.trt, source = "ec")
  super_pop <- pi_fun(pop_rct, pop_ec, mu.x2 = c(0, mu.x2.ec), sigma.x2 = c(1, sigma.x2.ec), lambda = lambda)

  if (pop_only == T) {
    # to learn about the method
    weighted_pop <- super_pop %>%
      mutate(
        Z = ifelse(source == "trial", 1, 0),
        pate.w1 = 1 / pi, pate.w0 = 1 / (1 - pi), patt.w1 = 1, patt.w0 = pi / (1 - pi), 
        patec.w1 = (1 - pi)/pi, patec.w0 = 1, pato.w1 = 1 - pi, pato.w0 = pi
      ) %>%
      mutate(
        pate.nw1 = pate.w1 / sum(Z * pate.w1), pate.nw0 = pate.w0 / sum((1 - Z) * pate.w0),
        patt.nw1 = patt.w1 / sum(Z * patt.w1), patt.nw0 = patt.w0 / sum((1 - Z) * patt.w0),
        patec.nw1 = patec.w1 / sum(Z * patec.w1), patec.nw0 = patec.w0 / sum((1 - Z) * patec.w0),
        pato.nw1 = pato.w1 / sum(Z * pato.w1), pato.nw0 = pato.w0 / sum((1 - Z) * pato.w0)
      ) %>%
      mutate(w.pate = Z * pate.nw1 + (1 - Z) * pate.nw0, w.patt = Z * patt.nw1 + (1 - Z) * patt.nw0, 
             w.patec = Z * patec.nw1 + (1 - Z) * patec.nw0, w.pato = Z * pato.nw1 + (1 - Z) * pato.nw0) %>%
      mutate(phi1 = phi1, phi2 = phi2, mu.x2.ec = mu.x2.ec, sigma.x2.ec = sigma.x2.ec, b.trt = b.trt, lambda = lambda)

    return(weighted_pop)
  } else {
    # for real simulation
    output <- NULL
    for (type in c("PATE", "PATT", "PATEC", "PATO")) {
      # corresponding to ATI, ATT, ATEC and ATO in manuscript
      output <- rbind(output, output_fun(super_pop, N1 = N1, alloc_ratio = alloc_ratio, N2 = N2, type = type, 
                                         phi1 = phi1, phi2 = phi2, mu.x2.ec = mu.x2.ec, sigma.x2.ec = sigma.x2.ec))
    }
    return(output)
  }
}
