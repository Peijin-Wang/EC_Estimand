# EC Estimand Paper -- case study
# analysis
#
# Peijin Wang
# 2025-6-9


library(tidyverse)
library(PSweight)
library(rms)
# library(forestplot)
library(ggforestplot)
library(patchwork) # merge plots
library(ggsci) #color
library(latex2exp)
library(kableExtra)
library(Hmisc) # weighted quantile


setwd("\\\\prd-dcri-smb03.azure.dhe.duke.edu/dcri/ct/pcori115jia/staff/pw134")
source("functions.R")

setwd("\\\\prd-dcri-smb03.azure.dhe.duke.edu/dcri/ct/pcori115jia/staff/bz91")
raw_data <- read.csv("PoolDS_De_idtf.csv")

# read data
data <- raw_data %>% mutate(A = ifelse(ITTARM %in% c(2,4), 1, 0), # treatment
                            Z = ifelse(ITTARM %in% c(1,2,4), 1, 0), # data source
                            phase = ifelse(ITTARM %in% c(1,2), 1, 2)) %>% # phase 1/2 of LIMIT-JIA
  mutate(BLNJNT_0 = ifelse(BLNJNT == 0, 1, 0),
         BLNJNT_1 = ifelse(BLNJNT == 1, 1, 0),
         BLNJNT_2 = ifelse(BLNJNT == 2, 1, 0),
         BLNJNT_3 = ifelse(BLNJNT == 3, 1, 0),
         FEMALE = ifelse(SEX == 2, 1, 0)) %>% 
 select(A, Z, phase, AGEREF, FEMALE, DIAGDYS, UVSTRAT, BLINJECT, BLNJNT,
        BLNJNT_0, BLNJNT_1, BLNJNT_2, BLNJNT_3, JADAS10, PHYSGLBL)


# phase II data
trial <- data %>% filter(phase == 2) %>% filter(Z == 1)
ec <- data %>% filter(phase == 2) %>% filter(Z == 0)

# Table 1 ####
covariates_binary <- c("BLNJNT_0", "BLNJNT_1", "BLNJNT_2", "BLNJNT_3","BLINJECT", "FEMALE", "UVSTRAT")
covariates_cont <- c("JADAS10", "PHYSGLBL")
covariates_count <- c("AGEREF", "DIAGDYS")
covariates <- c("BLNJNT_0", "BLNJNT_1", "BLNJNT_2", "BLNJNT_3","JADAS10",
                "PHYSGLBL","BLINJECT", 
                "AGEREF", "FEMALE", "DIAGDYS", "UVSTRAT")

table1 <- reduce(list(table1_fun(trial, ec, covariates, "raw", "raw"),
                      table1_fun(trial, ec, covariates, "PATE", "combine"),
                      table1_fun(trial, ec, covariates, "PATT", "combine"),
                      table1_fun(trial, ec, covariates, "PATC", "combine"),
                      table1_fun(trial, ec, covariates, "PATO", "combine")),
                 function(x,y) left_join(x,y, by = "covariate"))

kb_data <- table1 %>%
  mutate(
    sw_X_PATE = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_PATE * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_PATE, sw_X_q25_PATE, sw_X_q75_PATE),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_PATE, sw_X_q25_PATE, sw_X_q75_PATE),
      TRUE                             ~ NA_character_
    ),
    sw_X_PATT = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_PATT * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_PATT, sw_X_q25_PATT, sw_X_q75_PATT),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_PATT, sw_X_q25_PATT, sw_X_q75_PATT),
      TRUE                             ~ NA_character_
    ),
    sw_X_PATC = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_PATC * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_PATC, sw_X_q25_PATC, sw_X_q75_PATC),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_PATC, sw_X_q25_PATC, sw_X_q75_PATC),
      TRUE                             ~ NA_character_
    ),
    sw_X_PATO = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_PATO * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_PATO, sw_X_q25_PATO, sw_X_q75_PATO),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_PATO, sw_X_q25_PATO, sw_X_q75_PATO),
      TRUE                             ~ NA_character_
    ),
    trial = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", trial_mean * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 trial_q50, trial_q25, trial_q75),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 trial_q50, trial_q25, trial_q75),
      TRUE                             ~ NA_character_
    ),
    EC = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", EC_mean * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 EC_q50, EC_q25, EC_q75),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 EC_q50, EC_q25, EC_q75),
      TRUE                             ~ NA_character_
    )
  ) %>%
  select(covariate, trial, EC, sw_X_PATE, sw_X_PATT, sw_X_PATC, sw_X_PATO)

kable(kb_data, format = "latex", booktabs = T, digits = 2, linesep  = "")

# Table 1 by treatment ####

table1 <- reduce(list(table1_fun(trial, ec, covariates, "raw", "raw"),
                      table1_fun(trial, ec, covariates, "PATE", "by Z"),
                      table1_fun(trial, ec, covariates, "PATT", "by Z"),
                      table1_fun(trial, ec, covariates, "PATC", "by Z"),
                      table1_fun(trial, ec, covariates, "PATO", "by Z")),
                 function(x,y) left_join(x,y, by = "covariate"))

kb_data <- table1 %>%
  mutate(
    sw_X_0_PATE = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_0_PATE * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_0_PATE, sw_X_q25_0_PATE, sw_X_q75_0_PATE),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_0_PATE, sw_X_q25_0_PATE, sw_X_q75_0_PATE),
      TRUE                             ~ NA_character_
    ),
    sw_X_1_PATE = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_1_PATE * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_1_PATE, sw_X_q25_1_PATE, sw_X_q75_1_PATE),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_1_PATE, sw_X_q25_1_PATE, sw_X_q75_1_PATE),
      TRUE                             ~ NA_character_
    ),
    sw_X_0_PATT = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_0_PATT * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_0_PATT, sw_X_q25_0_PATT, sw_X_q75_0_PATT),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_0_PATT, sw_X_q25_0_PATT, sw_X_q75_0_PATT),
      TRUE                             ~ NA_character_
    ),
    sw_X_1_PATT = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_1_PATT * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_1_PATT, sw_X_q25_1_PATT, sw_X_q75_1_PATT),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_1_PATT, sw_X_q25_1_PATT, sw_X_q75_1_PATT),
      TRUE                             ~ NA_character_
    ),
    sw_X_0_PATC = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_0_PATC * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_0_PATC, sw_X_q25_0_PATC, sw_X_q75_0_PATC),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_0_PATC, sw_X_q25_0_PATC, sw_X_q75_0_PATC),
      TRUE                             ~ NA_character_
    ),
    sw_X_1_PATC = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_1_PATC * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_1_PATC, sw_X_q25_1_PATC, sw_X_q75_1_PATC),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_1_PATC, sw_X_q25_1_PATC, sw_X_q75_1_PATC),
      TRUE                             ~ NA_character_
    ),
    sw_X_0_PATO = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_0_PATO * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_0_PATO, sw_X_q25_0_PATO, sw_X_q75_0_PATO),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_0_PATO, sw_X_q25_0_PATO, sw_X_q75_0_PATO),
      TRUE                             ~ NA_character_
    ),
    sw_X_1_PATO = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", sw_X_mean_1_PATO * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 sw_X_q50_1_PATO, sw_X_q25_1_PATO, sw_X_q75_1_PATO),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 sw_X_q50_1_PATO, sw_X_q25_1_PATO, sw_X_q75_1_PATO),
      TRUE                             ~ NA_character_
    ),
    trial = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", trial_mean * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 trial_q50, trial_q25, trial_q75),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 trial_q50, trial_q25, trial_q75),
      TRUE                             ~ NA_character_
    ),
    EC = case_when(
      covariate %in% covariates_binary ~ sprintf("%.0f", EC_mean * 100),
      covariate %in% covariates_cont   ~ sprintf("%.2f (%.2f, %.2f)",
                                                 EC_q50, EC_q25, EC_q75),
      covariate %in% covariates_count  ~ sprintf("%.0f (%.0f, %.0f)",
                                                 EC_q50, EC_q25, EC_q75),
      TRUE                             ~ NA_character_
    )
  ) %>%
  select(covariate, trial, EC, sw_X_1_PATE, sw_X_0_PATE, sw_X_1_PATT, sw_X_0_PATT, 
         sw_X_1_PATC, sw_X_0_PATC, sw_X_1_PATO, sw_X_0_PATO)

kable(kb_data, format = "latex", booktabs = T, digits = 2, linesep  = "")










#### Backup Code ####
# Figure 2(b) mimic ####
# learn about distributions (not included in the paper)
# Original/Weighted Density Comparison -- stabilized weight

w_pate <- data_sw_fun(trial, ec, "PATE") %>% mutate(norm_sw = Z*sw1/sum(Z*w1) + (1-Z)*sw0/sum((1-Z)*w0)) %>% 
  mutate(type = "PATE")
w_patt <- data_sw_fun(trial, ec, "PATT") %>% mutate(norm_sw = Z*sw1/sum(Z*w1) + (1-Z)*sw0/sum((1-Z)*w0)) %>% 
  mutate(type = "PATT")
w_patc <- data_sw_fun(trial, ec, "PATC") %>% mutate(norm_sw = Z*sw1/sum(Z*w1) + (1-Z)*sw0/sum((1-Z)*w0)) %>% 
  mutate(type = "PATC")
w_pato <- data_sw_fun(trial, ec, "PATO") %>% mutate(norm_sw = Z*sw1/sum(Z*w1) + (1-Z)*sw0/sum((1-Z)*w0)) %>% 
  mutate(type = "PATO")

dt_plot_raw <- rbind(trial, ec)
dt_plot_weighted <- rbind(w_pate, w_patt, w_patc, w_pato)

## continuous variables
colors <- pal_locuszoom()(7)
cont_plot_fun <- function(covariate, x_lab) {
  p1 <- ggplot() + 
    geom_density(data = dt_plot_raw, aes(x = !!sym(covariate), color = "Pool"), linewidth = 1) +
    geom_density(data = dt_plot_raw %>% filter(Z == 1), aes(x = !!sym(covariate), color = "Trial"), linewidth = 1) +
    geom_density(data = dt_plot_raw %>% filter(Z == 0), aes(x = !!sym(covariate), color = "EC"), linewidth = 1) +
    scale_color_manual(values = c("Pool" = colors[1],"Trial" = colors[2],"EC" = colors[3]),
                      breaks = c("Pool", "Trial", "EC") ,name = NULL) +
    labs(y = "Density", x = x_lab) +
    theme_bw() +
    theme(legend.position = "bottom")
    
  
  p2 <- ggplot() +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATE"), aes(x = !!sym(covariate), weight = norm_sw, color = "ATI"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATT"), aes(x = !!sym(covariate), weight = norm_sw, color = "ATT"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATC"), aes(x = !!sym(covariate), weight = norm_sw, color = "ATEC"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATO"), aes(x = !!sym(covariate), weight = norm_sw, color = "ATO"),
                 linewidth = 1) +
    scale_color_manual(values = c("ATI" = colors[4], "ATT" = colors[5], "ATEC" = colors[6], "ATO" = colors[7]),
                       breaks = c("ATI", "ATT", "ATEC", "ATO"),
                       name = NULL) +
    labs(y = "Density", x = x_lab) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  p1 + p2
}



# continuous variables

pdf("Z:/ct/pcori115jia/staff/pw134/results/20250429_cont_plot.pdf", height = 12, width = 8)
(cont_plot_fun("AGEREF", "Age (AGEREF)")) /
  (cont_plot_fun("DIAGDYS", "Diagnosis Days (DIAGDYS)")) /
  (cont_plot_fun("JADAS10", "Disease Activity Score (JADAS10)")) /
  (cont_plot_fun("PHYSGLBL", "Physician Global Score (PHYSGLBL)"))
dev.off()

## other categorical variables
covariate <- "FEMALE"
level1 <- "Female"
level0 <- "Male"

cat_plot_fun <- function(covariate, levels, fill_lab) {
  level1 <- levels[1]
  level0 <- levels[2]
  
  dt_plot_weighted_sub <- dt_plot_weighted %>% group_by(type) %>% 
    summarise(X = sum(norm_sw * !!sym(covariate))) %>% 
    mutate(X0 = 1 - X) %>% 
    pivot_longer(X:X0, names_to = "covariate", values_to = "prop") %>% 
    mutate(level = ifelse(covariate == "X", level1, level0))
  
  dt_plot_raw_sub <- rbind(
    dt_plot_raw %>% group_by(Z) %>% summarise(X = mean(!!sym(covariate))) %>% 
      mutate(X0 = 1 - X, type = ifelse(Z == 1, "Trial", "EC")) %>% select(-Z),
    dt_plot_raw %>% summarise(X = mean(!!sym(covariate))) %>% 
      mutate(X0 = 1 - X, type = "Pool")) %>% 
    pivot_longer(X:X0,  names_to = "covariate", values_to = "prop") %>% 
    mutate(level = ifelse(covariate == "X", level1, level0))
  
  
  dt_plot_long <- rbind(dt_plot_weighted_sub, dt_plot_raw_sub)
  dt_plot_long$type <- factor(dt_plot_long$type, levels = c("Pool", "Trial", "EC", "PATE", "PATT", "PATC", "PATO"),
                              labels = c("Pool", "Trial", "EC", "ATI", "ATT", "ATEC", "ATO"))
  
  ggplot() +
    geom_col(data = dt_plot_long, aes(x = prop, y = type, fill = level), alpha = 0.7) +
    scale_fill_locuszoom() +
    labs(x = "Proportion", y = "", fill = fill_lab) +
    theme_bw() + theme(legend.position = "bottom")
}

# Physician global score

physglbl_plot_fun <- function(covariate = "PHYSGLBL", x_lab = "Physician Global Score (PHYSGLBL)") {
  p1 <- ggplot() + 
    geom_density(data = dt_plot_raw, aes(x = !!sym(covariate), color = "Pool"), 
                 linewidth = 1) +
    geom_density(data = dt_plot_raw %>% filter(Z == 1), aes(x = !!sym(covariate), color = "Trial"), linewidth = 1) +
    geom_density(data = dt_plot_raw %>% filter(Z == 0), aes(x = !!sym(covariate), color = "EC"), linewidth = 1) +
    scale_color_manual(values = c("Pool" = colors[1],"Trial" = colors[2],"EC" = colors[3]),
                       breaks = c("Pool", "Trial", "EC") ,name = NULL) +
    scale_x_continuous(limits = c(-5, 10)) +
    labs(y = "Density", x = x_lab) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  
  p2 <- ggplot() +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATE"), aes(x = !!sym(covariate), weight = norm_sw, color = "ATI"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATT"), aes(x = !!sym(covariate), weight = norm_sw, color = "ATT"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATC"), aes(x = !!sym(covariate), weight = norm_sw, color = "ATEC"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATO"), aes(x = !!sym(covariate), weight = norm_sw, color = "ATO"),
                 linewidth = 1) +
    scale_color_manual(values = c("ATI" = colors[4], "ATT" = colors[5], "ATEC" = colors[6], "ATO" = colors[7]),
                       breaks = c("ATI", "ATT", "ATEC", "ATO"),
                       name = NULL) +
    scale_x_continuous(limits = c(-5, 10)) +
    labs(y = "Density", x = x_lab) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  p1 + p2
}

pdf("\\\\prd-dcri-smb03.azure.dhe.duke.edu/dcri/ct/pcori115jia/staff/pw134/results/20250505_physglbl_plot.pdf", 
    height = 4, width = 8)
physglbl_plot_fun()
dev.off()

log_physglbl_plot_fun <- function(covariate = "PHYSGLBL", x_lab = "Log Physician Global Score (PHYSGLBL)") {
  p1 <- ggplot() + 
    geom_density(data = dt_plot_raw, aes(x = log(!!sym(covariate)), color = "Pool"), 
                 linewidth = 1) +
    geom_density(data = dt_plot_raw %>% filter(Z == 1), aes(x = log(!!sym(covariate)), color = "Trial"), linewidth = 1) +
    geom_density(data = dt_plot_raw %>% filter(Z == 0), aes(x = log(!!sym(covariate)), color = "EC"), linewidth = 1) +
    scale_color_manual(values = c("Pool" = colors[1],"Trial" = colors[2],"EC" = colors[3]),
                       breaks = c("Pool", "Trial", "EC") ,name = NULL) +
    scale_x_continuous(limits = c(-5, 5)) +
    labs(y = "Density", x = x_lab) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  
  p2 <- ggplot() +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATE"), aes(x = log(!!sym(covariate)), weight = norm_sw, color = "ATI"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATT"), aes(x = log(!!sym(covariate)), weight = norm_sw, color = "ATT"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATC"), aes(x = log(!!sym(covariate)), weight = norm_sw, color = "ATEC"),
                 linewidth = 1) +
    geom_density(data = dt_plot_weighted %>% filter(type == "PATO"), aes(x = log(!!sym(covariate)), weight = norm_sw, color = "ATO"),
                 linewidth = 1) +
    scale_color_manual(values = c("ATI" = colors[4], "ATT" = colors[5], "ATEC" = colors[6], "ATO" = colors[7]),
                       breaks = c("ATI", "ATT", "ATEC", "ATO"),
                       name = NULL) +
    scale_x_continuous(limits = c(-5, 5)) +
    labs(y = "Density", x = x_lab) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  p1 + p2
}

pdf("\\\\prd-dcri-smb03.azure.dhe.duke.edu/dcri/ct/pcori115jia/staff/pw134/results/20250505_log_physglbl_plot.pdf", 
    height = 4, width = 8)
log_physglbl_plot_fun()
dev.off()



## BLNJNT
dt_plot_weighted_sub <- dt_plot_weighted %>% group_by(type) %>% 
  summarise(X0 = sum(norm_sw * BLNJNT_0), X1 = sum(norm_sw * BLNJNT_1),
            X2 = sum(norm_sw * BLNJNT_2), X3 = sum(norm_sw * BLNJNT_3)) %>% 
  pivot_longer(X0:X3, names_to = "covariate", values_to = "prop") %>% 
  mutate(level = ifelse(covariate == "X0", 0, ifelse(covariate == "X1", 1,
                                                     ifelse(covariate == "X2", 2, 3))))

dt_plot_raw_sub <- rbind(
  dt_plot_raw %>% group_by(Z) %>% 
    summarise(X0 = mean(BLNJNT_0), X1 = mean(BLNJNT_1),
              X2 = mean(BLNJNT_2), X3 = mean(BLNJNT_3)) %>% 
    mutate(type = ifelse(Z == 1, "Trial", "EC")) %>% select(-Z),
  dt_plot_raw %>% summarise(X0 = mean(BLNJNT_0), X1 = mean(BLNJNT_1),
                            X2 = mean(BLNJNT_2), X3 = mean(BLNJNT_3)) %>% 
    mutate(type = "Pool")) %>% 
  pivot_longer(X0:X3,  names_to = "covariate", values_to = "prop") %>% 
  mutate(level = ifelse(covariate == "X0", 0, ifelse(covariate == "X1", 1,
                                                     ifelse(covariate == "X2", 2, 3))))


dt_plot_long <- rbind(dt_plot_weighted_sub, dt_plot_raw_sub)
dt_plot_long$type <- factor(dt_plot_long$type, levels = c("Pool", "Trial", "EC", "PATE", "PATT", "PATC", "PATO"),
                            labels = c("Pool", "Trial", "EC", "ATI", "ATT", "ATEC", "ATO"))

p_BLNJNT <- ggplot() +
  geom_col(data = dt_plot_long, aes(x = prop, y = type, fill = as.factor(level)), alpha = 0.7) +
  scale_fill_locuszoom() +
  labs(x = "Proportion", y = "", fill = "# Active Joints (BLNJNT)") +
  theme_bw() + theme(legend.position = "bottom")





pdf("Z:/ct/pcori115jia/staff/pw134/results/20250429_cat_plot.pdf", height = 12, width = 8)
cat_plot_fun("FEMALE", c("Female", "Male"), "Sex") /
  cat_plot_fun("UVSTRAT", c("Uveitis Risk Higher", "Uveitis Risk Lower"), "Uveitis Strata (UVSTRAT)") /
  cat_plot_fun("BLINJECT", c("Yes", "No"), "Glucocorticoids (BLINJECT)") /
  p_BLNJNT
dev.off()

# Number Active Joints bar plot
dt_plot_long <- rbind(dt_plot_weighted_sub, dt_plot_raw_sub) %>% filter(type %in% c("PATE", "PATT", "PATC", "PATO"))
dt_plot_long$type <- factor(dt_plot_long$type, levels = c("PATC", "PATE", "PATO", "PATT"),
                            labels = c("ATI", "ATT", "ATEC", "ATO"))
dt_plot_long$level <- factor(dt_plot_long$level, levels = 0:3)

pdf("\\\\prd-dcri-smb03.azure.dhe.duke.edu/dcri/ct/pcori115jia/staff/pw134/results/20250505_blnjnt_plot.pdf", 
    height = 4, width = 8)
ggplot() +
  geom_col(data = dt_plot_long, aes(x = prop, y = type, fill = level),position = position_stack(reverse = TRUE), alpha = 0.7) +
  scale_fill_locuszoom() +
  labs(x = "Proportion", y = "", fill = "# Active Joints (BLNJNT)") +
  theme_bw() + theme(legend.position = "bottom")
dev.off()


