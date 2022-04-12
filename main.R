library(magrittr)
library(simulator)
library(tidyverse)
library(optimx)

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

###############################################################
### TYPE 1 ERROR SIMULATIONS WITH CORRELATED DATA AND E-M #####
###############################################################
my_cores <- parallel::detectCores()

my_nsim = 200
new_simulation(name = "vary_n_type1",
              label = "varying by n type1 simulations") %>%
  generate_model(make_data_correlated_v2, seed = 123,
                 n = as.list(c(100,50,30)), 
                 beta0 = 0.1, 
                 beta1 = 0,
                 Mi_mu = 2.4, 
                 Mi_var = 0.5, 
                 epsilon = 0,
                 gamma0 = -2,
                 gamma1 = 0.1, 
                 vary_along = "n") %>%
  simulate_from_model(nsim = my_nsim/my_cores, index = 1:my_cores) %>% 
  run_method(c(GLM_Rao,david_method2), parallel = list(socket_names = my_cores)) -> type1_123_methods

type1_123_methods %>% evaluate(list(pvalue)) -> type1_123_evals

type1_123_evals %>% evals %>% as.data.frame -> vary_123
type1_123_evals %>% model %>% as.data.frame -> model_123
type1_123_dataframe <- dplyr::right_join(model_123,  vary_123, by = c("name" = "Model"))
type1_123_dataframe$p05 <- dplyr::if_else(type1_2022_dataframe$pvalue < 0.05, 1, 0) 




