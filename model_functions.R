## @knitr models
require(magrittr)
require(tidyverse)
library(mgcv)
require(mgcv)
make_data <- function(n, 
                      Mi_mu,
                      Mi_var, 
                      beta0, 
                      beta1,
                      g1_mu,
                      g1_var,
                      epsilon){
  make_table <- function(n, 
                         Mi_mu,
                         Mi_var, 
                         beta0,
                         beta1,
                         g1_mu,
                         g1_var,
                         epsilon){
    X1 <- rbinom(n,1,0.5)
    Mi_log <- rnorm(n,Mi_mu,Mi_var) #total number of reads
    Mi <- exp(Mi_log)
    coverage = Mi
    log_coverage = log(coverage)
    gamma1 <- g1_mu
    beta.table <- cbind(X1,beta0,beta1) %>% as_tibble
    beta.table$lambda_props = exp(beta.table$beta0 + beta.table$beta1*beta.table$X1)/(1+exp(beta.table$beta0+beta.table$beta1*beta.table$X1))
    beta.table$lambda <- rbinom(n, 1, beta.table$lambda_props)
    gamma.table <- cbind(beta.table,coverage,log_coverage,gamma1) %>% as_tibble
    gamma.table$Y_prop <- ifelse(gamma.table$lambda == 1, 1-exp(-gamma1*gamma.table$coverage), epsilon)
    gamma.table$Y <- rbinom(n, 1, gamma.table$Y_prop)
    return(gamma.table)
  }
  new_model(name = "generate_data", 
            label = sprintf("Parameters (n = %s, Mi_mu = %s, Mi_var = %s, beta0 = %s, beta1 = %s, g1_mu = %s, g1_var = %s, epsilon = %s)", 
                            n, 
                            Mi_mu,
                            Mi_var, 
                            beta0,
                            beta1,
                            g1_mu,
                            g1_var,
                            epsilon), 
            params = list(n = n,
                          Mi_mu = Mi_mu,
                          Mi_var = Mi_var, 
                          beta0 = beta0,
                          beta1 = beta1, 
                          g1_mu = g1_mu,
                          g1_var = g1_var,
                          epsilon = epsilon), 
            simulate = function(n, 
                                Mi_mu,
                                Mi_var, 
                                beta0,
                                beta1,
                                g1_mu,
                                g1_var,
                                epsilon,
                                nsim){
              sim_list <- replicate(nsim, make_table(n, 
                                                     Mi_mu,
                                                     Mi_var, 
                                                     beta0,
                                                     beta1,
                                                     g1_mu,
                                                     g1_var,
                                                     epsilon), simplify = F)})
}



# Yay I think I finished up this data generation
make_data_correlated <- function(n, 
                                 beta0, 
                                 beta1,
                                 Mi_mu, 
                                 Mi_var,
                                 gamma,
                                 epsilon){
  
  make_table <- function(n, 
                         beta0, 
                         beta1,
                         Mi_mu, 
                         Mi_var, 
                         gamma,
                         epsilon){
    Mi_log <- rnorm(n, Mi_mu, Mi_var) #total number of reads
    Mi <- exp(Mi_log)
    coverage = Mi
    log_coverage = log(coverage)
    my_coverage_data <- cbind(log_coverage,coverage) %>% as_tibble()
    
    my_coverage_data %>% 
      mutate(X1 = rbinom(n(),1,boot::inv.logit(gamma*coverage)), 
             prob_X = boot::inv.logit(gamma*coverage)) -> my_coverage_data
    
    # When gamma = 0 then prob_X = 0.5 
    # Let's try cases where prob_X = 0.6, 0.7, 0.8
    
    beta.table <- cbind(my_coverage_data,beta0,beta1) %>% as_tibble()
    
    # gamma1 <- g1_mu
    beta.table$lambda_props = exp(beta.table$beta0 + beta.table$beta1*beta.table$X1)/(1+exp(beta.table$beta0+beta.table$beta1*beta.table$X1))
    beta.table$lambda <- rbinom(n, 1, beta.table$lambda_props)
    
    gamma.table <- beta.table %>% as_tibble
    
    ###### Read in the Arimizu data 
    
    readRDS("logit_DRR102664.RDS") -> logit_data
    gam_model <- gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")
    # pred <- predict(gam_model, type = "response", newdata = gamma.table) # output the probabilities 
    gamma.table$Y_prop <- ifelse(gamma.table$lambda == 1, predict(gam_model, type = "response", newdata = gamma.table), epsilon)
    gamma.table$Y <- rbinom(n, 1, gamma.table$Y_prop)
    return(gamma.table)
  }
  
  new_model(name = "generate_data", 
            label = sprintf("Parameters (n = %s, beta0 = %s, beta1 = %s,Mi_mu = %s, Mi_var = %s, gamma = %s, epsilon = %s)", 
                            n, 
                            beta0, 
                            beta1,
                            Mi_mu, 
                            Mi_var,
                            gamma,
                            epsilon), 
            params = list(n = n,
                          beta0 = beta0,
                          beta1 = beta1,
                          Mi_mu = Mi_mu, 
                          Mi_var = Mi_mu,
                          gamma = gamma,
                          epsilon = epsilon), 
            simulate = function(n, 
                                beta0, 
                                beta1, 
                                Mi_mu, 
                                Mi_var,
                                gamma,
                                epsilon,
                                nsim){
              sim_list <- replicate(nsim, make_table(n, 
                                                     beta0, 
                                                     beta1,
                                                     Mi_mu, 
                                                     Mi_var, 
                                                     gamma,
                                                     epsilon), simplify = F)})
}

make_data_correlated_v2 <- function(n, 
                                    Mi_mu,
                                    Mi_var, 
                                    beta0, 
                                    beta1,
                                    gamma0,
                                    gamma1,
                                    epsilon){
  
  make_table <- function(n, 
                         Mi_mu,
                         Mi_var, 
                         beta0, 
                         beta1,
                         gamma0,
                         gamma1,
                         epsilon){
    
    Mi_log <- rnorm(n, Mi_mu, Mi_var) #total number of reads
    Mi <- exp(Mi_log)
    coverage = Mi
    log_coverage = Mi_log
    my_coverage_data <- cbind(log_coverage,coverage) %>% as_tibble()
    
    my_coverage_data %>% 
      mutate(numerator = exp(gamma0+gamma1*coverage), 
             denominator = 1+ exp(gamma0 + gamma1*coverage), 
             prob_x = numerator/denominator,
             prob_X = boot::inv.logit(gamma0 + gamma1*coverage), 
             X1 = rbinom(n,1,prob_X)) -> my_coverage_data
    
    # ggplot(my_coverage_data, aes(x = coverage, y = X1)) + geom_point()
    # ggplot(my_coverage_data, aes(x = coverage, y = prob_X)) + geom_point()
    
    # I want x to be negative for lower M_i's but more positive for higher M_i's 
    
    # When gamma = 0 then prob_X = 0.5 
    # Let's try cases where prob_X = 0.6, 0.7, 0.8
    
    beta.table <- cbind(my_coverage_data,beta0,beta1) %>% as_tibble()
    
    beta.table %>% 
      mutate(lambda_props = exp(beta.table$beta0 + beta.table$beta1*beta.table$X1)/(1+exp(beta.table$beta0+beta.table$beta1*beta.table$X1))) -> beta.table
    
    beta.table$lambda <- rbinom(n, 1, beta.table$lambda_props)
    
    gamma.table <- beta.table %>% as_tibble
    
    ###### Read in the Arimizu data 
    
    readRDS("logit_DRR102664.RDS") -> logit_data
    
    gam_model <- gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")
    # pred <- predict(gam_model, type = "response", newdata = gamma.table) # output the probabilities 
    gamma.table %>% 
      mutate(Y_prop = ifelse(gamma.table$lambda == 1, predict(gam_model, type = "response", newdata = gamma.table), epsilon)) -> gamma.table
    gamma.table$Y <- rbinom(n,1,gamma.table$Y_prop)
    #ggplot(gamma.table, aes(x = coverage, y = Y, color = as.factor(X1))) + geom_point()
    return(gamma.table)
    
  }
  
  new_model(name = "generate_data", 
            label = sprintf("Parameters (n = %s, Mi_mu = %s, Mi_var = %s, beta0 = %s, beta1 = %s, gamma0 = %s, gamma1 = %s, epsilon = %s)", 
                            n, 
                            Mi_mu,
                            Mi_var, 
                            beta0, 
                            beta1,
                            gamma0,
                            gamma1,
                            epsilon), 
            params = list(n = n,
                          Mi_mu = Mi_mu,
                          Mi_var = Mi_var, 
                          beta0 = beta0,
                          beta1 = beta1,
                          gamma0 = gamma0,
                          gamma1 = gamma1,
                          epsilon = epsilon), 
            simulate = function(n, 
                                Mi_mu,
                                Mi_var, 
                                beta0, 
                                beta1,
                                gamma0, 
                                gamma1,
                                epsilon,
                                nsim){
              sim_list <- replicate(nsim, make_table(n, 
                                                     Mi_mu,
                                                     Mi_var, 
                                                     beta0, 
                                                     beta1,
                                                     gamma0,
                                                     gamma1,
                                                     epsilon), simplify = F)})
}


make_data_section4 <- function(n, 
                               Mi_mu,
                               Mi_var, 
                               beta0, 
                               beta1,
                               epsilon){
  
  make_table <- function(n, 
                         Mi_mu,
                         Mi_var, 
                         beta0, 
                         beta1,
                         epsilon){
    X1 <- rbinom(n,1,0.5)
    Mi_log <- rnorm(n, Mi_mu, Mi_var) #total number of reads
    Mi <- exp(Mi_log)
    coverage = Mi
    log_coverage = log(coverage)
    # gamma1 <- g1_mu
    beta.table <- cbind(X1,beta0,beta1) %>% as_tibble
    beta.table$lambda_props = exp(beta.table$beta0 + beta.table$beta1*beta.table$X1)/(1+exp(beta.table$beta0+beta.table$beta1*beta.table$X1))
    beta.table$lambda <- rbinom(n, 1, beta.table$lambda_props)
    
    gamma.table <- cbind(beta.table,coverage,log_coverage) %>% as_tibble
    
    ###### Read in the Arimizu data 
    
    readRDS("logit_DRR102664.RDS") -> logit_data
    gam_model <- gam(presence ~ s(coverage), family = binomial, data = logit_data, method = "REML")
    # pred <- predict(gam_model, type = "response", newdata = gamma.table) # output the probabilities 
    gamma.table$Y_prop <- ifelse(gamma.table$lambda == 1, predict(gam_model, type = "response", newdata = gamma.table), epsilon)
    gamma.table$Y <- rbinom(n, 1, gamma.table$Y_prop)
    return(gamma.table)
  }
  
  new_model(name = "generate_data", 
            label = sprintf("Parameters (n = %s, Mi_mu = %s, Mi_var = %s, beta0 = %s, beta1 = %s, epsilon = %s)", 
                            n, 
                            Mi_mu,
                            Mi_var, 
                            beta0, 
                            beta1, 
                            epsilon), 
            params = list(n = n,
                          Mi_mu = Mi_mu,
                          Mi_var = Mi_var, 
                          beta0 = beta0,
                          beta1 = beta1, 
                          epsilon = epsilon), 
            simulate = function(n, 
                                Mi_mu,
                                Mi_var, 
                                beta0, 
                                beta1, 
                                epsilon,
                                nsim){
              sim_list <- replicate(nsim, make_table(n, 
                                                     Mi_mu,
                                                     Mi_var, 
                                                     beta0, 
                                                     beta1, 
                                                     epsilon), simplify = F)})
}
