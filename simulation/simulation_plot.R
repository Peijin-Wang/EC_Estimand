# EC Estimand Project -- Simulation Results Output
#
# Peijin Wang
# 2025-6-9

library(tidyverse)
library(latex2exp)
library(ggsci)
library(patchwork)
library(kableExtra)

rm(list = ls())
setwd("/Users/wangpeijin/Desktop/PhD Research/EC Estimand Project/1 Simulation/20250512/results")
files <- list.files(pattern = "\\.csv$")
data_list <- lapply(files, read.csv)
out <- bind_rows(data_list)
out$type <- as.factor(out$type)

setwd("/Users/wangpeijin/Desktop/PhD Research/EC Estimand Project/1 Simulation/20250512/sim_output")

# Simulation Summary Data ####
summary <- out %>% group_by(type, N1, N11, N10, N2, lambda, phi, mu.x2.ec, sigma.x2.ec) %>% 
  summarise(bias = mean(estimator-estimand), mse = mean((estimator-estimand)^2),
            m.estimand = mean(estimand), m.estimator = mean(estimator)) %>% ungroup() %>% 
  mutate(EC = case_when(mu.x2.ec == 0 & sigma.x2.ec == 1 ~ "EC1",
                        mu.x2.ec == 0.5 & sigma.x2.ec == 1 ~ "EC2",
                        mu.x2.ec == 1 & sigma.x2.ec == 1 ~ "EC3",
                        mu.x2.ec == 2 & sigma.x2.ec == 1 ~ "EC4",
                        mu.x2.ec == 0 & sigma.x2.ec != 1 ~ "EC5",
                        mu.x2.ec == 0.5 & sigma.x2.ec != 1 ~ "EC6",
                        mu.x2.ec == 1 & sigma.x2.ec != 1 ~ "EC7",
                        mu.x2.ec == 2 & sigma.x2.ec != 1 ~ "EC8",
                        TRUE ~ "Mistake") )


# Bias and MSE Plots ####
plot_fun <- function(n10=100, ec="EC1", y_var="bias", ylim = c(0,1), title){
  plot_dt <- summary %>% mutate(phi = factor(phi), type = factor(type),
                                ratio = ifelse(N10 != 0, paste0("N10:N2=",1,":", N2/N10), 
                                               paste0("N1:N2=",1,":", N2/N1))) %>% 
    filter(N10 == n10 & EC == ec)
  
  plot_dt$type <- factor(plot_dt$type, levels = c("PATE","PATT","PATEC","PATO"))
  levels(plot_dt$type) <- c("ATI","ATT","ATEC","ATO")
  plot_dt$ratio <- factor(plot_dt$ratio, levels = c("N10:N2=1:1", "N10:N2=1:3", "N10:N2=1:10",
                                                    "N1:N2=1:1", "N1:N2=1:3", "N1:N2=1:10"))
  if(n10==100){
    plot_dt$setting <- case_when(plot_dt$ratio == "N10:N2=1:1" ~ "setting 1-3",
                                 plot_dt$ratio == "N10:N2=1:3" ~ "setting 4-6",
                                 TRUE ~ "setting 7-9")
  }else if(n10==50){
    plot_dt$setting <- case_when(plot_dt$ratio == "N10:N2=1:1" ~ "setting 10-12",
                                 plot_dt$ratio == "N10:N2=1:3" ~ "setting 13-15",
                                 TRUE ~ "setting 16-18")
  }else if(n10==0){
    plot_dt$setting <- case_when(plot_dt$ratio == "N1:N2=1:1" ~ "setting 19-21",
                                 plot_dt$ratio == "N1:N2=1:3" ~ "setting 22-24",
                                 TRUE ~ "setting 25-27")
  }
  
  if(y_var == "bias"){
    p <- ggplot(data=plot_dt, aes(x = phi, y = .data[[y_var]], group = type, color = type)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.8) +
      geom_point() + geom_line() +
      facet_grid(.~setting+ratio) +
      scale_color_locuszoom() +
      labs(x = TeX("HTE ($phi_1=phi_2$)"), y = "Bias", color = NULL, title = title) +
      ylim(ylim) +
      theme_bw() +
      theme(plot.title = element_text(size = 20),        # Title font size
            plot.subtitle = element_text(size = 18),     # Title font size
            axis.title.x = element_text(size = 16),      # X-axis title font size
            axis.title.y = element_text(size = 16),      # Y-axis title font size
            axis.text.x = element_text(size = 14),       # X-axis text font size
            axis.text.y = element_text(size = 14),       # Y-axis text font size
            legend.title = element_text(size = 16),      # Legend title font size
            legend.text = element_text(size = 14),       # Legend text font size
            strip.text.x = element_text(size = 14)
      )
  }
  
  if(y_var == "mse"){
    p <- ggplot(data=plot_dt, aes(x = phi, y = .data[[y_var]], group = type, color = type)) +
      geom_point() + geom_line() +
      facet_grid(.~setting+ratio) +
      scale_color_locuszoom() +
      labs(x = TeX("HTE ($phi_1=phi_2$)"), y = "MSE", color = NULL, title = title) +
      ylim(ylim) +
      theme_bw() +
      theme(plot.title = element_text(size = 20),        # Title font size
            plot.subtitle = element_text(size = 18),        # Title font size
            axis.title.x = element_text(size = 16),      # X-axis title font size
            axis.title.y = element_text(size = 16),      # Y-axis title font size
            axis.text.x = element_text(size = 14),       # X-axis text font size
            axis.text.y = element_text(size = 14),       # Y-axis text font size
            legend.title = element_text(size = 16),      # Legend title font size
            legend.text = element_text(size = 14),        # Legend text font size
            strip.text.x = element_text(size = 14)
      )
  }
  return(p)
}

plots <- list(plot_fun(n10=100, ec = "EC1", y_var = "bias", ylim = c(-0.1,0.1), title = "EC1 with X2 ~ N(0,1)"),
              plot_fun(n10=100, ec = "EC2", y_var = "bias", ylim = c(-0.1,0.1), title = "EC2 with X2 ~ N(0.5,1)"),
              plot_fun(n10=100, ec = "EC3", y_var = "bias", ylim = c(-0.5,0.5), title = "EC3 with X2 ~ N(1,1)"),
              plot_fun(n10=100, ec = "EC4", y_var = "bias", ylim = c(-0.5,0.5), title = "EC4 with X2 ~ N(2,1)"),
              plot_fun(n10=100, ec = "EC5", y_var = "bias", ylim = c(-0.3,0.2), title = "EC5 with X2 ~ N(0,1.5)"),
              plot_fun(n10=100, ec = "EC6", y_var = "bias", ylim = c(-0.3,0.2), title = "EC6 with X2 ~ N(0.5,1.5)"),
              plot_fun(n10=100, ec = "EC7", y_var = "bias", ylim = c(-0.8,0.3), title = "EC7 with X2 ~ N(1,1.5)"),
              plot_fun(n10=100, ec = "EC8", y_var = "bias", ylim = c(-0.8,0.3), title = "EC8 with X2 ~ N(2,1.5)"))


pdf(file = "bias1.pdf", width = 12, height = 14)
wrap_plots(plots, ncol = 2) +  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
dev.off()

plots <- list(plot_fun(n10=50, ec = "EC1", y_var = "bias", ylim = c(-0.1,0.1), title = "EC1 with X2 ~ N(0,1)"),
              plot_fun(n10=50, ec = "EC2", y_var = "bias", ylim = c(-0.1,0.1), title = "EC2 with X2 ~ N(0.5,1)"),
              plot_fun(n10=50, ec = "EC3", y_var = "bias", ylim = c(-0.5,0.5), title = "EC3 with X2 ~ N(1,1)"),
              plot_fun(n10=50, ec = "EC4", y_var = "bias", ylim = c(-0.5,0.5), title = "EC4 with X2 ~ N(2,1)"),
              plot_fun(n10=50, ec = "EC5", y_var = "bias", ylim = c(-0.2,0.2), title = "EC5 with X2 ~ N(0,1.5)"),
              plot_fun(n10=50, ec = "EC6", y_var = "bias", ylim = c(-0.2,0.2), title = "EC6 with X2 ~ N(0.5,1.5)"),
              plot_fun(n10=50, ec = "EC7", y_var = "bias", ylim = c(-0.8,0.5), title = "EC7 with X2 ~ N(1,1.5)"),
              plot_fun(n10=50, ec = "EC8", y_var = "bias", ylim = c(-0.8,0.5), title = "EC8 with X2 ~ N(2,1.5)"))


pdf(file = "bias2.pdf", width = 12, height = 14)
wrap_plots(plots, ncol = 2) +  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
dev.off()


plots <- list(plot_fun(n10=0, ec = "EC1", y_var = "bias", ylim = c(-0.1,0.1), title = "EC1 with X2 ~ N(0,1)"),
              plot_fun(n10=0, ec = "EC2", y_var = "bias", ylim = c(-0.1,0.1), title = "EC2 with X2 ~ N(0.5,1)"),
              plot_fun(n10=0, ec = "EC3", y_var = "bias", ylim = c(-0.5,0.5), title = "EC3 with X2 ~ N(1,1)"),
              plot_fun(n10=0, ec = "EC4", y_var = "bias", ylim = c(-0.5,0.5), title = "EC4 with X2 ~ N(2,1)"),
              plot_fun(n10=0, ec = "EC5", y_var = "bias", ylim = c(-0.3,0.2), title = "EC5 with X2 ~ N(0,1.5)"),
              plot_fun(n10=0, ec = "EC6", y_var = "bias", ylim = c(-0.3,0.2), title = "EC6 with X2 ~ N(0.5,1.5)"),
              plot_fun(n10=0, ec = "EC7", y_var = "bias", ylim = c(-0.8,0.5), title = "EC7 with X2 ~ N(1,1.5)"),
              plot_fun(n10=0, ec = "EC8", y_var = "bias", ylim = c(-0.8,0.5), title = "EC8 with X2 ~ N(2,1.5)"))


pdf(file = "bias3.pdf", width = 12, height = 14)
wrap_plots(plots, ncol = 2) +  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
dev.off()



## plot for ICHPS poster
plots <- list(plot_fun(n11=150, ec = "EC1", y_var = "bias", ylim = c(-0.5,0.5), title = "EC1 with X2 ~ N(0,1)"),
              plot_fun(n11=150, ec = "EC4", y_var = "bias", ylim = c(-0.5,0.5), title = "EC4 with X2 ~ N(2,1)"),
              plot_fun(n11=150, ec = "EC6", y_var = "bias", ylim = c(-0.5,0.5), title = "EC6 with X2 ~ N(0.5,1.5)"),
              plot_fun(n11=150, ec = "EC8", y_var = "bias", ylim = c(-0.5,0.5), title = "EC8 with X2 ~ N(2,1.5)"))


pdf(file = "bias_poster.pdf", width = 16, height = 8)
wrap_plots(plots, ncol = 2) +  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
dev.off()

# MSE
plots <- list(plot_fun(n10=100, ec = "EC1", y_var = "mse",  ylim = c(0,0.1), title = "EC1 with X2 ~ N(0,1)"),
              plot_fun(n10=100, ec = "EC2", y_var = "mse",  ylim = c(0,0.1), title = "EC2 with X2 ~ N(0.5,1)"),
              plot_fun(n10=100, ec = "EC3", y_var = "mse",  ylim = c(0,0.8), title = "EC3 with X2 ~ N(1,1)"),
              plot_fun(n10=100, ec = "EC4", y_var = "mse",  ylim = c(0,0.8), title = "EC4 with X2 ~ N(2,1)"),
              plot_fun(n10=100, ec = "EC5", y_var = "mse",  ylim = c(0,0.1), title = "EC5 with X2 ~ N(0,1.5)"),
              plot_fun(n10=100, ec = "EC6", y_var = "mse",  ylim = c(0,0.1), title = "EC6 with X2 ~ N(0.5,1.5)"),
              plot_fun(n10=100, ec = "EC7", y_var = "mse",  ylim = c(0,1.2), title = "EC7 with X2 ~ N(1,1.5)"),
              plot_fun(n10=100, ec = "EC8", y_var = "mse",  ylim = c(0,1.2), title = "EC8 with X2 ~ N(2,1.5)"))


pdf(file = "mse1.pdf", width = 12, height = 14)
wrap_plots(plots, ncol = 2) +  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
dev.off()

plots <- list(plot_fun(n10=50, ec = "EC1", y_var = "mse", ylim = c(0,0.1), title = "EC1 with X2 ~ N(0,1)"),
              plot_fun(n10=50, ec = "EC2", y_var = "mse", ylim = c(0,0.1), title = "EC2 with X2 ~ N(0.5,1)"),
              plot_fun(n10=50, ec = "EC3", y_var = "mse", ylim = c(0,0.8), title = "EC3 with X2 ~ N(1,1)"),
              plot_fun(n10=50, ec = "EC4", y_var = "mse", ylim = c(0,0.8), title = "EC4 with X2 ~ N(2,1)"),
              plot_fun(n10=50, ec = "EC5", y_var = "mse", ylim = c(0,0.1), title = "EC5 with X2 ~ N(0,1.5)"),
              plot_fun(n10=50, ec = "EC6", y_var = "mse", ylim = c(0,0.1), title = "EC6 with X2 ~ N(0.5,1.5)"),
              plot_fun(n10=50, ec = "EC7", y_var = "mse", ylim = c(0,0.8), title = "EC7 with X2 ~ N(1,1.5)"),
              plot_fun(n10=50, ec = "EC8", y_var = "mse", ylim = c(0,0.8), title = "EC8 with X2 ~ N(2,1.5)"))


pdf(file = "mse2.pdf", width = 12, height = 14)
wrap_plots(plots, ncol = 2) +  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
dev.off()

plots <- list(plot_fun(n10=0, ec = "EC1", y_var = "mse", ylim = c(0,0.1), title = "EC1 with X2 ~ N(0,1)"),
              plot_fun(n10=0, ec = "EC2", y_var = "mse", ylim = c(0,0.1), title = "EC2 with X2 ~ N(0.5,1)"),
              plot_fun(n10=0, ec = "EC3", y_var = "mse", ylim = c(0,0.8), title = "EC3 with X2 ~ N(1,1)"),
              plot_fun(n10=0, ec = "EC4", y_var = "mse", ylim = c(0,0.8), title = "EC4 with X2 ~ N(2,1)"),
              plot_fun(n10=0, ec = "EC5", y_var = "mse", ylim = c(0,0.1), title = "EC5 with X2 ~ N(0,1.5)"),
              plot_fun(n10=0, ec = "EC6", y_var = "mse", ylim = c(0,0.1), title = "EC6 with X2 ~ N(0.5,1.5)"),
              plot_fun(n10=0, ec = "EC7", y_var = "mse", ylim = c(0,1.0), title = "EC7 with X2 ~ N(1,1.5)"),
              plot_fun(n10=0, ec = "EC8", y_var = "mse", ylim = c(0,1.0), title = "EC8 with X2 ~ N(2,1.5)"))


pdf(file = "mse3.pdf", width = 12, height = 14)
wrap_plots(plots, ncol = 2) +  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
dev.off()



# Estimand Table ####
tb.estimand <- summary %>% select(type, phi, lambda, m.estimand, EC) %>% pivot_wider(names_from = type, values_from = m.estimand)
tb.estimand <- tb.estimand %>% arrange(lambda, phi, EC) %>% select(phi, EC, lambda, PATE, PATT, PATEC, PATO)
kable(tb.estimand %>% filter(lambda==round(2/3,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.estimand %>% filter(lambda==round(2/5,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.estimand %>% filter(lambda==round(1/6,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.estimand %>% filter(lambda==round(4/5,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.estimand %>% filter(lambda==round(4/7,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.estimand %>% filter(lambda==round(2/7,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.estimand %>% filter(lambda==round(1/2,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.estimand %>% filter(lambda==round(1/4,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.estimand %>% filter(lambda==round(1/11,16)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "",digits = 2)

# Bias Table ####
tb.bias <- summary %>% select(type, phi, lambda, bias, EC) %>% pivot_wider(names_from = type, values_from = bias)
tb.bias <- tb.bias %>% arrange(lambda, phi, EC) %>% select(phi, EC, lambda, PATE, PATT, PATEC, PATO)
kable(tb.bias %>% filter(lambda==round(2/3,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.bias %>% filter(lambda==round(2/5,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.bias %>% filter(lambda==round(1/6,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.bias %>% filter(lambda==round(4/5,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.bias %>% filter(lambda==round(4/7,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.bias %>% filter(lambda==round(2/7,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.bias %>% filter(lambda==round(1/2,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.bias %>% filter(lambda==round(1/4,15)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)
kable(tb.bias %>% filter(lambda==round(1/11,16)) %>% select(-phi,-lambda), format = "latex", booktabs = TRUE, linesep = "", digits = 2)


