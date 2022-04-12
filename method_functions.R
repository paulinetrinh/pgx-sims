## @knitr methods
require(optimx)
require(isotone)

LRT_noepsilon_clust <- new_method("LRT_noepsilon_clust", "LRT_noepsilon_clust", 
                                  method = function(model,draw){
                                    restricted_function_noepsilon = function(par, my_data){
                                      e = 0 # e is parameter argument 1 
                                      vBeta = c(par[1],0)
                                      vGamma = par[2]
                                      Y0 <- my_data %>% dplyr::filter(Y == 0)
                                      Y1 <- my_data %>% dplyr::filter(Y == 1)
                                      
                                      Mi0 = as.matrix(Y0['coverage'])
                                      #Mi0 = log(Mi0)
                                      # add an intercept to the predictor variables
                                      #Mi0 = cbind(1, Mi0)
                                      mX0 = as.matrix(Y0['X1'])
                                      mX0 = cbind(1, mX0)
                                      
                                      Mi1 = as.matrix(Y1['coverage'])
                                      # Mi1 = log(Mi1) # modeling the log coverage
                                      # Mi1 = cbind(1, Mi1)
                                      mX1 = as.matrix(Y1['X1'])
                                      mX1 = cbind(1, mX1)
                                      
                                      LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                      (1-(1-exp(-(Mi0 %*% vGamma))))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                        sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                  (1/(1+exp(-(mX1 %*% vBeta))))*(1-exp(-(Mi1 %*% vGamma)))))
                                      return(LL)
                                    }
                                    unrestricted_function_noepsilon = function(par, my_data){
                                      e = 0 # e is parameter argument 1 
                                      vBeta = c(par[1],par[2]) # my beta vector is parameter arguments 2 and 3 
                                      vGamma = par[3] # my gamma vector is parameter arguments 4
                                      Y0 <- my_data %>% dplyr::filter(Y == 0)
                                      Y1 <- my_data %>% dplyr::filter(Y == 1)
                                      
                                      Mi0 = as.matrix(Y0['coverage']) # I'm going to model the regular coverage
                                      # Mi0 = log(Mi0)
                                      # add an intercept to the predictor variables
                                      # Mi0 = cbind(1, Mi0)
                                      mX0 = as.matrix(Y0['X1'])
                                      mX0 = cbind(1, mX0)
                                      
                                      Mi1 = as.matrix(Y1['coverage'])
                                      # Mi1 = log(Mi1) # modeling the log coverage
                                      # Mi1 = cbind(1, Mi1)
                                      mX1 = as.matrix(Y1['X1'])
                                      mX1 = cbind(1, mX1)
                                      
                                      LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                      (1-(1-exp(-(Mi0*vGamma))))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                        sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                  ((1-exp(-(Mi1*vGamma)))/(1+exp(-(mX1 %*% vBeta))))))
                                      return(LL)
                                    }
                                    # Make a for loop that plugs in a set of parameters and then 
                                    # plugs in those to optim 
                                    # then returns the parameter and likelihood
                                    
                                    Y0 <- draw %>% dplyr::filter(Y == 0)
                                    Y1 <- draw %>% dplyr::filter(Y == 1)
                                    Mi0 = as.matrix(Y0['coverage'])
                                    #Mi0 = log(Mi0)
                                    Mi0_mean = mean(Mi0)
                                    
                                    Mi1 = as.matrix(Y1['coverage'])
                                    #Mi1 = log(Mi1)
                                    Mi1_mean = mean(Mi1) 
                                    
                                    gamma_init <- (-1/Mi1_mean)*log((Mi0_mean)/(Mi0_mean+Mi1_mean))
                                    # gamma_par = 0.11
                                    
                                    # take the maximum log-likelihood row out of the tibble
                                    checkWarning <- base::tryCatch(mylogit <- stats::glm(Y ~ X1, data = draw, family = "binomial"), 
                                                                   warning = function(w){
                                                                     if(base::grepl("algorithm did not converge",as.character(w)) == "TRUE"){
                                                                       z <- 1}
                                                                   }
                                    )
                                    
                                    if (length(checkWarning[[1]]) == 2){
                                      init_values <- c(mylogit$coefficients[[1]],mylogit$coefficients[[2]],gamma_init[[1]]) 
                                      init_values2 <- c(mylogit$coefficients[[1]],gamma_init[[1]])
                                      
                                      init_list <- list()
                                      init_list[[1]] <- init_values
                                      # init_list[[2]] <- c(0.0001,mylogit$coefficients[[1]]/2,mylogit$coefficients[[2]]/2,gamma_init[[1]]/2)
                                      init_list[[2]] <- c(mylogit$coefficients[[1]]/3,mylogit$coefficients[[2]]/3,gamma_init[[1]]/3)
                                      init_list[[3]] <- c(mylogit$coefficients[[1]]*2,mylogit$coefficients[[2]]*2,gamma_init[[1]]*2)
                                      # init_list[[5]] <- c(0.01,mylogit$coefficients[[1]]/2,mylogit$coefficients[[2]]/2,gamma_init[[1]]/2)
                                      init_list[[4]] <- c(0.1,0.1,0.1)
                                      init_list[[5]] <- c(0.00001,0.00001,0.00001)
                                      init_list[[6]] <- c(0.0001,0.0001,0.0001)
                                      init_list[[7]] <- c(0,0,1)
                                      init_list[[8]] <- c(0,0,15)
                                      init_list[[9]] <- c(0,0,100)
                                      init_list[[10]] <- c(0,0,200)
                                      init_list[[11]] <- c(0,0,300)
                                      init_list[[12]] <- c(0,0,500)
                                      
                                      init_list2 <- list()
                                      init_list2[[1]] <- init_values2
                                      init_list2[[2]] <- c(mylogit$coefficients[[1]]/3,gamma_init[[1]]/3)
                                      init_list2[[3]] <- c(mylogit$coefficients[[1]]*2,gamma_init[[1]]*2)
                                      # init_list2[[5]] <- c(0.001,mylogit$coefficients[[1]]/2,gamma_init[[1]]/2)
                                      init_list2[[4]] <- c(0.1,0.1)
                                      init_list2[[5]] <- c(0.00001,0.00001)
                                      init_list2[[6]] <- c(0.0001,0.0001)
                                      init_list2[[7]] <- c(0,1)
                                      init_list2[[8]] <- c(0,15)
                                      init_list2[[9]] <- c(0,100)
                                      init_list2[[10]] <- c(0,200)
                                      init_list2[[11]] <- c(0,300)
                                      init_list2[[12]] <- c(0,500)
                                      
                                    } else {
                                      init_list <- list() 
                                      init_list[[1]] <- c(0.1,0.1,0.1)
                                      init_list[[2]] <- c(0.00001,0.00001,0.00001)
                                      init_list[[3]] <- c(0.0001,0.0001,0.0001)
                                      init_list[[4]] <- c(0.1,0.1,gamma_init[[1]])
                                      init_list[[5]] <- c(1,1,gamma_init[[1]]*2)
                                      init_list[[6]] <- c(0,0,1)
                                      init_list[[7]] <- c(0,0,15)
                                      init_list[[8]] <- c(0,0,100)
                                      init_list[[9]] <- c(0,0,200)
                                      init_list[[10]] <- c(0,0,300)
                                      init_list[[11]] <- c(0,0,500)
                                      
                                      init_list2 <- list()
                                      init_list2[[1]] <- c(0.1,0.1)
                                      init_list2[[2]] <- c(0.00001,0.00001)
                                      init_list2[[3]] <- c(0.0001,0.0001)
                                      init_list2[[4]] <- c(0.1,gamma_init[[1]])
                                      init_list2[[5]] <- c(1,gamma_init[[1]]*2)
                                      init_list2[[6]] <- c(0,1)
                                      init_list2[[7]] <- c(0,15)
                                      init_list2[[8]] <- c(0,100)
                                      init_list2[[9]] <- c(0,200)
                                      init_list2[[10]] <- c(0,300)
                                      init_list2[[11]] <- c(0,500)
                                    }
                                    
                                    restricted_list <- list()
                                    for (i in 1:length(init_list2)){  
                                      restricted_model = optimx::optimx(init_list2[[i]], restricted_function_noepsilon, my_data = draw, 
                                                                        method = 'L-BFGS-B',
                                                                        lower=c(-Inf,0),
                                                                        upper=c(Inf, Inf),
                                                                        control = list(trace = 0, fnscale= -1)) 
                                      cbind(restricted_model$p1, restricted_model$p2,restricted_model$value) -> restricted_list[[i]]
                                    }  
                                    do.call(rbind, restricted_list) %>% 
                                      dplyr::as_tibble() %>% 
                                      dplyr::filter(!is.na(V1)) %>% 
                                      dplyr::distinct() -> restricted_tibble
                                    restricted_tibble[which.max(restricted_tibble$V3),] -> restricted_results 
                                    if (dim(restricted_results)[[1]] == 0){
                                      print("restricted model no epsilon didn't converge")
                                      restricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA)
                                      LRT = NA
                                      pvalue = NA
                                    } else {
                                      print("restricted model no epsilon completed!")
                                    }
                                    unrestricted_list <- list() 
                                    for (j in 1:length(init_list)){
                                      unrestricted_model = optimx::optimx(init_list[[j]], unrestricted_function_noepsilon, my_data = draw, 
                                                                          method = 'L-BFGS-B',
                                                                          lower= c(-Inf,-Inf, 0),
                                                                          upper=c(Inf, Inf, Inf),
                                                                          control = list(trace = 0, fnscale= -1))
                                      cbind(unrestricted_model$p1, unrestricted_model$p2,unrestricted_model$p3,unrestricted_model$value) -> unrestricted_list[[j]]
                                    }
                                    do.call(rbind, unrestricted_list) %>% 
                                      dplyr::as_tibble() %>% 
                                      dplyr::filter(!is.na(V1)) %>% 
                                      dplyr::distinct() -> unrestricted_tibble
                                    unrestricted_tibble[which.max(unrestricted_tibble$V4),] -> unrestricted_results 
                                    if (dim(unrestricted_results)[[1]] == 0){
                                      print("unrestricted model no epsilon didn't converge")
                                      unrestricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA)
                                      LRT = NA 
                                      pvalue = NA
                                    } else {
                                      print("unrestricted model no epsilon completed!")
                                      LRT = 2*abs(unrestricted_results$V4 - restricted_results$V3)
                                      LRT_old = 2*(unrestricted_results$V4 - restricted_results$V3)
                                      pvalue <- stats::pchisq(LRT, df = 1, lower.tail = FALSE)
                                    }
                                    
                                    output <- list("beta0" = unrestricted_results$V1,
                                                   "beta1" = unrestricted_results$V2,
                                                   "gamma1" = unrestricted_results$V3,
                                                   "pvalue" = pvalue, 
                                                   "LL_unrestricted" = unrestricted_results$V4, 
                                                   "LL_restricted" = restricted_results$V3, 
                                                   "LRT_abs" = LRT, 
                                                   "LRT_noabs" = LRT_old)
                                  }) 
# Create a quality function for 1-e^-(M_i*gamma)
quality_function <- function(quality_vector, gamma_vector){
  quality_function_value <- 1-exp(-(quality_vector %*% gamma_vector))
  return(quality_function_value)
}

quality_function_logit <- function(quality_vector, gamma_vector){
  quality_function_value <- exp(quality_vector %*% gamma_vector)/(1+exp(quality_vector %*% gamma_vector))
  return(quality_function_value)
}
# Creating a log likelihood function to use: 
calculate_LL <- function(input_e, input_vBeta, input_vGamma, input_mX0, input_Mi0, input_mX1, input_Mi1){
  LL <- sum(log((1-input_e)*exp(-(input_mX0 %*% input_vBeta))/(1+exp(-(input_mX0 %*% input_vBeta))) + # when Y = 0 
                  (1-(quality_function(input_Mi0, input_vGamma)))*(1/(1+exp(-(input_mX0 %*% input_vBeta)))))) + 
    sum(log(input_e*exp(-(input_mX1 %*% input_vBeta))/(1+exp(-(input_mX1 %*% input_vBeta))) +  # when Y = 1
              (1/(1+exp(-(input_mX1 %*% input_vBeta))))*(quality_function(input_Mi1,input_vGamma))))
  return(LL)
}


LRT_noepsilon_clust_logit <- new_method("LRT_noepsilon_logit", "LRT_noepsilon_logit", 
                                        method = function(model,draw){
                                          restricted_function_noepsilon_logit = function(par, my_data){
                                            e = 0 # e is parameter argument 1 
                                            vBeta = c(par[1],0)
                                            vGamma = c(par[2],par[3])
                                            Y0 <- my_data %>% dplyr::filter(Y == 0)
                                            Y1 <- my_data %>% dplyr::filter(Y == 1)
                                            
                                            Mi0 = as.matrix(Y0['coverage'])
                                            #Mi0 = log(Mi0)
                                            # add an intercept to the predictor variables
                                            Mi0 = cbind(1, Mi0)
                                            mX0 = as.matrix(Y0['X1'])
                                            mX0 = cbind(1, mX0)
                                            
                                            Mi1 = as.matrix(Y1['coverage'])
                                            # Mi1 = log(Mi1) # modeling the log coverage
                                            Mi1 = cbind(1, Mi1)
                                            mX1 = as.matrix(Y1['X1'])
                                            mX1 = cbind(1, mX1)
                                            
                                            LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                            (1-(quality_function_logit(quality_vector = Mi0, gamma_vector = vGamma)))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                              sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                        (1/(1+exp(-(mX1 %*% vBeta))))*(quality_function_logit(quality_vector = Mi1, gamma_vector = vGamma))))
                                            return(LL)
                                          }
                                          unrestricted_function_noepsilon_logit = function(par, my_data){
                                            e = 0 # e is parameter argument 1 
                                            vBeta = c(par[1],par[2]) # my beta vector is parameter arguments 2 and 3 
                                            vGamma = c(par[3],par[4]) # my gamma vector is parameter arguments 4
                                            Y0 <- my_data %>% dplyr::filter(Y == 0)
                                            Y1 <- my_data %>% dplyr::filter(Y == 1)
                                            
                                            Mi0 = as.matrix(Y0['coverage']) # I'm going to model the regular coverage
                                            # Mi0 = log(Mi0)
                                            # add an intercept to the predictor variables
                                            Mi0 = cbind(1, Mi0)
                                            mX0 = as.matrix(Y0['X1'])
                                            mX0 = cbind(1, mX0)
                                            
                                            Mi1 = as.matrix(Y1['coverage'])
                                            # Mi1 = log(Mi1) # modeling the log coverage
                                            Mi1 = cbind(1, Mi1)
                                            mX1 = as.matrix(Y1['X1'])
                                            mX1 = cbind(1, mX1)
                                            
                                            LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                            (1-(quality_function_logit(quality_vector = Mi0, gamma_vector = vGamma)))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                              sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                        (1/(1+exp(-(mX1 %*% vBeta))))*(quality_function_logit(quality_vector = Mi1, gamma_vector = vGamma))))
                                            return(LL)
                                          }
                                          # Make a for loop that plugs in a set of parameters and then 
                                          # plugs in those to optim 
                                          # then returns the parameter and likelihood
                                          
                                          #Y0 <- draw %>% dplyr::filter(Y == 0)
                                          #Y1 <- draw %>% dplyr::filter(Y == 1)
                                          #Mi0 = as.matrix(Y0['coverage'])
                                          #Mi0 = log(Mi0)
                                          #Mi0_mean = mean(Mi0)
                                          
                                          #Mi1 = as.matrix(Y1['coverage'])
                                          #Mi1 = log(Mi1)
                                          #Mi1_mean = mean(Mi1) 
                                          
                                          # Mi_data <- draw %>% dplyr::filter(lambda == 1)
                                          # What should my initializing points be for gamma? 
                                          
                                          checkWarning <- base::tryCatch(my_gammalogit <- stats::glm(Y ~ coverage, data = draw, family = "binomial"), 
                                                                         warning = function(w){
                                                                           if(base::grepl("algorithm did not converge",as.character(w)) == "TRUE"){
                                                                             z <- 1}
                                                                         }
                                          )
                                          #gamma_init <- (-1/Mi1_mean)*log((Mi0_mean)/(Mi0_mean+Mi1_mean))
                                          # gamma_par = 0.11
                                          
                                          # take the maximum log-likelihood row out of the tibble
                                          checkWarning <- base::tryCatch(mylogit <- stats::glm(Y ~ X1, data = draw, family = "binomial"), 
                                                                         warning = function(w){
                                                                           if(base::grepl("algorithm did not converge",as.character(w)) == "TRUE"){
                                                                             z <- 1}
                                                                         }
                                          )
                                          
                                          if (length(checkWarning[[1]]) == 2){
                                            init_values <- c(mylogit$coefficients[[1]],mylogit$coefficients[[2]],my_gammalogit$coefficients[[1]],my_gammalogit$coefficients[[2]]) 
                                            init_values2 <- c(mylogit$coefficients[[1]],my_gammalogit$coefficients[[1]],my_gammalogit$coefficients[[2]])
                                            
                                            init_list <- list()
                                            init_list[[1]] <- init_values
                                            # init_list[[2]] <- c(0.0001,mylogit$coefficients[[1]]/2,mylogit$coefficients[[2]]/2,gamma_init[[1]]/2)
                                            # init_list[[2]] <- c(mylogit$coefficients[[1]]/3,mylogit$coefficients[[2]]/3,gamma_init[[1]]/3)
                                            # init_list[[3]] <- c(mylogit$coefficients[[1]]*2,mylogit$coefficients[[2]]*2,gamma_init[[1]]*2)
                                            # init_list[[5]] <- c(0.01,mylogit$coefficients[[1]]/2,mylogit$coefficients[[2]]/2,gamma_init[[1]]/2)
                                            init_list[[2]] <- c(0.1,0.1,0.1,0.1)
                                            init_list[[3]] <- c(0.00001,0.00001,0.00001,0.00001)
                                            init_list[[4]] <- c(0.0001,0.0001,0.0001,0.00001)
                                            init_list[[5]] <- c(1,1,1,1)
                                            #init_list[[6]] <- c(10,10,10,10)
                                            # init_list[[7]] <- c(10,10,100)
                                            
                                            init_list2 <- list()
                                            init_list2[[1]] <- init_values2
                                            #init_list2[[2]] <- c(mylogit$coefficients[[1]]/3,gamma_init[[1]]/3)
                                            #init_list2[[3]] <- c(mylogit$coefficients[[1]]*2,gamma_init[[1]]*2)
                                            # init_list2[[5]] <- c(0.001,mylogit$coefficients[[1]]/2,gamma_init[[1]]/2)
                                            init_list2[[2]] <- c(0.1,0.1,0.1)
                                            init_list2[[3]] <- c(0.00001,0.00001,0.00001)
                                            init_list2[[4]] <- c(0.0001,0.0001,0.0001)
                                            init_list2[[5]] <- c(1,1,1)
                                            #init_list2[[6]] <- c(10,10,10)
                                            # init_list2[[7]] <- c(10,100)
                                            
                                          } else {
                                            init_list <- list() 
                                            init_list[[1]] <- c(0.1,0.1,0.1,0.1)
                                            init_list[[2]] <- c(0.00001,0.00001,0.00001,0.00001)
                                            init_list[[3]] <- c(0.0001,0.0001,0.0001,0.0001)
                                            init_list[[4]] <- c(0.1,0.1,my_gammalogit$coefficients[[1]],my_gammalogit$coefficients[[2]])
                                            init_list[[5]] <- c(1,1,my_gammalogit$coefficients[[1]]*2,my_gammalogit$coefficients[[2]]*2)
                                            init_list[[6]] <- c(1,1,1,1)
                                            #init_list[[7]] <- c(10,10,10,10)
                                            #init_list[[8]] <- c(10,10,100,100)
                                            #init_list[[9]] <- c(10,10,200)
                                            #init_list[[10]] <- c(10,10,300)
                                            #init_list[[11]] <- c(10,10,500)
                                            
                                            init_list2 <- list()
                                            init_list2[[1]] <- c(0.1,0.1,0.1)
                                            init_list2[[2]] <- c(0.00001,0.00001,0.00001)
                                            init_list2[[3]] <- c(0.0001,0.0001,0.0001)
                                            init_list2[[4]] <- c(0.1,my_gammalogit$coefficients[[1]],my_gammalogit$coefficients[[2]])
                                            init_list2[[5]] <- c(1,my_gammalogit$coefficients[[1]]*2,my_gammalogit$coefficients[[2]]*2)
                                            init_list2[[6]] <- c(10,10,10)
                                            #init_list2[[7]] <- c(1,1,1)
                                            #init_list2[[8]] <- c(10,100,100)
                                            #init_list2[[9]] <- c(10,200)
                                            #init_list2[[10]] <- c(10,300)
                                            #init_list2[[11]] <- c(10,500)
                                          }
                                          restricted_list <- list()
                                          for (i in 1:length(init_list2)){  
                                            restricted_model = optimx::optimx(init_list2[[i]], restricted_function_noepsilon_logit, my_data = draw, 
                                                                              method = 'BFGS',
                                                                              control = list(trace = 0, fnscale= -1)) 
                                            cbind(restricted_model$p1, restricted_model$p2,restricted_model$p3, restricted_model$value) -> restricted_list[[i]]
                                          }  
                                          do.call(rbind, restricted_list) %>% 
                                            dplyr::as_tibble() %>% 
                                            dplyr::filter(!is.na(V1)) %>% 
                                            dplyr::distinct() -> restricted_tibble
                                          restricted_tibble[which.max(restricted_tibble$V4),] -> restricted_results 
                                          if (dim(restricted_results)[[1]] == 0){
                                            print("restricted model no epsilon didn't converge")
                                            restricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA)
                                            LRT = NA
                                            pvalue = NA
                                          } else {
                                            print("restricted model no epsilon completed!")
                                          }
                                          
                                          unrestricted_list <- list() 
                                          for (j in 1:length(init_list)){
                                            unrestricted_model = optimx::optimx(init_list[[j]], unrestricted_function_noepsilon_logit, my_data = draw, 
                                                                                method = 'BFGS',
                                                                                control = list(trace = 0, fnscale= -1))
                                            cbind(unrestricted_model$p1, unrestricted_model$p2,unrestricted_model$p3,unrestricted_model$p4,unrestricted_model$value) -> unrestricted_list[[j]]
                                          }
                                          do.call(rbind, unrestricted_list) %>% 
                                            dplyr::as_tibble() %>% 
                                            dplyr::filter(!is.na(V1)) %>% 
                                            dplyr::distinct() -> unrestricted_tibble
                                          unrestricted_tibble[which.max(unrestricted_tibble$V5),] -> unrestricted_results 
                                          if (dim(unrestricted_results)[[1]] == 0){
                                            print("unrestricted model no epsilon didn't converge")
                                            unrestricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                            LRT = NA 
                                            pvalue = NA
                                          } else {
                                            print("unrestricted model no epsilon completed!")
                                            LRT = 2*abs(unrestricted_results$V4 - restricted_results$V3)
                                            LRT_old = 2*(unrestricted_results$V4 - restricted_results$V3)
                                            pvalue <- stats::pchisq(LRT, df = 1, lower.tail = FALSE)
                                          }
                                          
                                          output <- list("beta0" = unrestricted_results$V1,
                                                         "beta1" = unrestricted_results$V2,
                                                         "gamma0" = unrestricted_results$V3,
                                                         "gamma1" = unrestricted_results$V4,
                                                         "pvalue" = pvalue, 
                                                         "LL_unrestricted" = unrestricted_results$V5, 
                                                         "LL_restricted" = restricted_results$V4, 
                                                         "LRT_abs" = LRT, 
                                                         "LRT_noabs" = LRT_old)
                                        }) 

LRT_noepsilon_clust_logit2 <- new_method("LRT_noepsilon_logit2", "LRT_noepsilon_logit2", 
                                         method = function(model,draw){
                                           restricted_function_noepsilon_logit = function(par, my_data){
                                             e = 0 # e is parameter argument 1 
                                             vBeta = c(par[1],0)
                                             vGamma = c(par[2],par[3])
                                             Y0 <- my_data %>% dplyr::filter(Y == 0)
                                             Y1 <- my_data %>% dplyr::filter(Y == 1)
                                             
                                             Mi0 = as.matrix(Y0['log_coverage'])
                                             #Mi0 = log(Mi0_1)
                                             # add an intercept to the predictor variables
                                             Mi0 = cbind(1, Mi0)
                                             mX0 = as.matrix(Y0['X1'])
                                             mX0 = cbind(1, mX0)
                                             
                                             Mi1 = as.matrix(Y1['log_coverage'])
                                             #Mi1 = log(Mi1_1) # modeling the log coverage
                                             Mi1 = cbind(1, Mi1)
                                             mX1 = as.matrix(Y1['X1'])
                                             mX1 = cbind(1, mX1)
                                             
                                             LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                             (1-(quality_function_logit(quality_vector = Mi0, gamma_vector = vGamma)))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                               sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                         (1/(1+exp(-(mX1 %*% vBeta))))*(quality_function_logit(quality_vector = Mi1, gamma_vector = vGamma))))
                                             return(LL)
                                           }
                                           unrestricted_function_noepsilon_logit = function(par, my_data){
                                             e = 0 # e is parameter argument 1 
                                             vBeta = c(par[1],par[2]) # my beta vector is parameter arguments 2 and 3 
                                             vGamma = c(par[3],par[4]) # my gamma vector is parameter arguments 4
                                             Y0 <- my_data %>% dplyr::filter(Y == 0)
                                             Y1 <- my_data %>% dplyr::filter(Y == 1)
                                             
                                             Mi0 = as.matrix(Y0['log_coverage']) # I'm going to model the regular coverage
                                             # Mi0 = log(Mi0_1)
                                             # add an intercept to the predictor variables
                                             Mi0 = cbind(1, Mi0)
                                             mX0 = as.matrix(Y0['X1'])
                                             mX0 = cbind(1, mX0)
                                             
                                             Mi1 = as.matrix(Y1['log_coverage'])
                                             #Mi1 = log(Mi1_1) # modeling the log coverage
                                             Mi1 = cbind(1, Mi1)
                                             mX1 = as.matrix(Y1['X1'])
                                             mX1 = cbind(1, mX1)
                                             
                                             LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                             (1-(quality_function_logit(quality_vector = Mi0, gamma_vector = vGamma)))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                               sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                         (1/(1+exp(-(mX1 %*% vBeta))))*(quality_function_logit(quality_vector = Mi1, gamma_vector = vGamma))))
                                             return(LL)
                                           }
                                           # Make a for loop that plugs in a set of parameters and then 
                                           # plugs in those to optim 
                                           # then returns the parameter and likelihood
                                           
                                           #Y0 <- draw %>% dplyr::filter(Y == 0)
                                           #Y1 <- draw %>% dplyr::filter(Y == 1)
                                           #Mi0 = as.matrix(Y0['coverage'])
                                           #Mi0 = log(Mi0)
                                           #Mi0_mean = mean(Mi0)
                                           
                                           #Mi1 = as.matrix(Y1['coverage'])
                                           #Mi1 = log(Mi1)
                                           #Mi1_mean = mean(Mi1) 
                                           
                                           # Mi_data <- draw %>% dplyr::filter(lambda == 1)
                                           # What should my initializing points be for gamma? 
                                           
                                           checkWarning <- base::tryCatch(my_gammalogit <- stats::glm(Y ~ log_coverage, data = draw, family = "binomial"), 
                                                                          warning = function(w){
                                                                            if(base::grepl("algorithm did not converge",as.character(w)) == "TRUE"){
                                                                              z <- 1}
                                                                          }
                                           )
                                           #gamma_init <- (-1/Mi1_mean)*log((Mi0_mean)/(Mi0_mean+Mi1_mean))
                                           # gamma_par = 0.11
                                           
                                           # take the maximum log-likelihood row out of the tibble
                                           checkWarning <- base::tryCatch(mylogit <- stats::glm(Y ~ X1, data = draw, family = "binomial"), 
                                                                          warning = function(w){
                                                                            if(base::grepl("algorithm did not converge",as.character(w)) == "TRUE"){
                                                                              z <- 1}
                                                                          }
                                           )
                                           
                                           if (length(checkWarning[[1]]) == 2){
                                             init_values <- c(mylogit$coefficients[[1]],mylogit$coefficients[[2]],my_gammalogit$coefficients[[1]],my_gammalogit$coefficients[[2]]) 
                                             init_values2 <- c(mylogit$coefficients[[1]],my_gammalogit$coefficients[[1]],my_gammalogit$coefficients[[2]])
                                             
                                             init_list <- list()
                                             init_list[[1]] <- init_values
                                             # init_list[[2]] <- c(0.0001,mylogit$coefficients[[1]]/2,mylogit$coefficients[[2]]/2,gamma_init[[1]]/2)
                                             # init_list[[2]] <- c(mylogit$coefficients[[1]]/3,mylogit$coefficients[[2]]/3,gamma_init[[1]]/3)
                                             # init_list[[3]] <- c(mylogit$coefficients[[1]]*2,mylogit$coefficients[[2]]*2,gamma_init[[1]]*2)
                                             # init_list[[5]] <- c(0.01,mylogit$coefficients[[1]]/2,mylogit$coefficients[[2]]/2,gamma_init[[1]]/2)
                                             init_list[[2]] <- c(0.1,0.1,0.1,0.1)
                                             init_list[[3]] <- c(0.00001,0.00001,0.00001,0.00001)
                                             init_list[[4]] <- c(0.0001,0.0001,0.0001,0.00001)
                                             init_list[[5]] <- c(1,1,1,1)
                                             #init_list[[6]] <- c(10,10,10,10)
                                             # init_list[[7]] <- c(10,10,100)
                                             
                                             init_list2 <- list()
                                             init_list2[[1]] <- init_values2
                                             #init_list2[[2]] <- c(mylogit$coefficients[[1]]/3,gamma_init[[1]]/3)
                                             #init_list2[[3]] <- c(mylogit$coefficients[[1]]*2,gamma_init[[1]]*2)
                                             # init_list2[[5]] <- c(0.001,mylogit$coefficients[[1]]/2,gamma_init[[1]]/2)
                                             init_list2[[2]] <- c(0.1,0.1,0.1)
                                             init_list2[[3]] <- c(0.00001,0.00001,0.00001)
                                             init_list2[[4]] <- c(0.0001,0.0001,0.0001)
                                             init_list2[[5]] <- c(1,1,1)
                                             #init_list2[[6]] <- c(10,10,10)
                                             # init_list2[[7]] <- c(10,100)
                                             
                                           } else {
                                             init_list <- list() 
                                             init_list[[1]] <- c(0.1,0.1,0.1,0.1)
                                             init_list[[2]] <- c(0.00001,0.00001,0.00001,0.00001)
                                             init_list[[3]] <- c(0.0001,0.0001,0.0001,0.0001)
                                             init_list[[4]] <- c(0.1,0.1,my_gammalogit$coefficients[[1]],my_gammalogit$coefficients[[2]])
                                             init_list[[5]] <- c(1,1,my_gammalogit$coefficients[[1]]*2,my_gammalogit$coefficients[[2]]*2)
                                             init_list[[6]] <- c(1,1,1,1)
                                             #init_list[[7]] <- c(10,10,10,10)
                                             #init_list[[8]] <- c(10,10,100,100)
                                             #init_list[[9]] <- c(10,10,200)
                                             #init_list[[10]] <- c(10,10,300)
                                             #init_list[[11]] <- c(10,10,500)
                                             
                                             init_list2 <- list()
                                             init_list2[[1]] <- c(0.1,0.1,0.1)
                                             init_list2[[2]] <- c(0.00001,0.00001,0.00001)
                                             init_list2[[3]] <- c(0.0001,0.0001,0.0001)
                                             init_list2[[4]] <- c(0.1,my_gammalogit$coefficients[[1]],my_gammalogit$coefficients[[2]])
                                             init_list2[[5]] <- c(1,my_gammalogit$coefficients[[1]]*2,my_gammalogit$coefficients[[2]]*2)
                                             init_list2[[6]] <- c(10,10,10)
                                             #init_list2[[7]] <- c(1,1,1)
                                             #init_list2[[8]] <- c(10,100,100)
                                             #init_list2[[9]] <- c(10,200)
                                             #init_list2[[10]] <- c(10,300)
                                             #init_list2[[11]] <- c(10,500)
                                           }
                                           restricted_list <- list()
                                           for (i in 1:length(init_list2)){  
                                             restricted_model = optimx::optimx(init_list2[[i]], restricted_function_noepsilon_logit, my_data = draw, 
                                                                               method = 'BFGS',
                                                                               control = list(trace = 0, fnscale= -1)) 
                                             cbind(restricted_model$p1, restricted_model$p2,restricted_model$p3, restricted_model$value) -> restricted_list[[i]]
                                           }  
                                           do.call(rbind, restricted_list) %>% 
                                             dplyr::as_tibble() %>% 
                                             dplyr::filter(!is.na(V1)) %>% 
                                             dplyr::distinct() -> restricted_tibble
                                           restricted_tibble[which.max(restricted_tibble$V4),] -> restricted_results 
                                           if (dim(restricted_results)[[1]] == 0){
                                             print("restricted model no epsilon didn't converge")
                                             restricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA)
                                             LRT = NA
                                             pvalue = NA
                                           } else {
                                             print("restricted model no epsilon completed!")
                                           }
                                           
                                           unrestricted_list <- list() 
                                           for (j in 1:length(init_list)){
                                             unrestricted_model = optimx::optimx(init_list[[j]], unrestricted_function_noepsilon_logit, my_data = draw, 
                                                                                 method = 'BFGS',
                                                                                 control = list(trace = 0, fnscale= -1))
                                             cbind(unrestricted_model$p1, unrestricted_model$p2,unrestricted_model$p3,unrestricted_model$p4,unrestricted_model$value) -> unrestricted_list[[j]]
                                           }
                                           do.call(rbind, unrestricted_list) %>% 
                                             dplyr::as_tibble() %>% 
                                             dplyr::filter(!is.na(V1)) %>% 
                                             dplyr::distinct() -> unrestricted_tibble
                                           unrestricted_tibble[which.max(unrestricted_tibble$V5),] -> unrestricted_results 
                                           if (dim(unrestricted_results)[[1]] == 0){
                                             print("unrestricted model no epsilon didn't converge")
                                             unrestricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                             LRT = NA 
                                             pvalue = NA
                                           } else {
                                             print("unrestricted model no epsilon completed!")
                                             LRT = 2*abs(unrestricted_results$V4 - restricted_results$V3)
                                             LRT_old = 2*(unrestricted_results$V4 - restricted_results$V3)
                                             pvalue <- stats::pchisq(LRT, df = 1, lower.tail = FALSE)
                                           }
                                           
                                           output <- list("beta0" = unrestricted_results$V1,
                                                          "beta1" = unrestricted_results$V2,
                                                          "gamma0" = unrestricted_results$V3,
                                                          "gamma1" = unrestricted_results$V4,
                                                          "pvalue" = pvalue, 
                                                          "LL_unrestricted" = unrestricted_results$V5, 
                                                          "LL_restricted" = restricted_results$V4, 
                                                          "LRT_abs" = LRT, 
                                                          "LRT_noabs" = LRT_old)
                                         }) 

## old method
david_method <- new_method("david_method","david_method",
                           method = function(model,draw){
                             ### Functions for the 3 hierarchical model components 
                             # Probability of Y = 1 given Lambda = 1
                             library(isotone)
                             Pr_Y1_Lambda1 <- function(my_fMi) {
                               my_Pr_Y1_Lambda1 <- exp(my_fMi)/(1+exp(my_fMi))
                               return(my_Pr_Y1_Lambda1)
                             }
                             # Probability of Lambda = 1 
                             Pr_Lambda1 <- function(my_data) {
                               mX1 = as.matrix(my_data['X1'])
                               mX1 = cbind(1, mX1)
                               my_data %>% 
                                 dplyr::select(beta0,beta1) %>% 
                                 dplyr::distinct() %>% 
                                 as.numeric() -> my_beta
                               my_Pr_Lambda1 <- exp(mX1 %*% my_beta)/(1+exp(mX1 %*% my_beta))
                               return(my_Pr_Lambda1)
                             }
                             # Probability of Y = 1 given Lambda = 0 
                             Pr_Y1_Lambda0 <- function(my_epsilon) {
                               #my_Pr_Y1_Lambda0 <- exp(my_epsilon)/(1+exp(my_epsilon))
                               # for now I'm going to set this parameterization to 0 
                               my_Pr_Y1_Lambda0 <- 0 
                               return(my_Pr_Y1_Lambda0)
                             }
                             
                             ### Functions for calculating the denominators of p_i 
                             # Prob of Y = 1 given parameters theta 
                             Pr_Y1_theta <- function(my_data){
                               my_Pr_Y1_theta <- Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data) + Pr_Y1_Lambda0(my_data$epsilon)*(1-Pr_Lambda1(my_data))
                               return(my_Pr_Y1_theta)
                             }
                             #Prob of Y = 0 given parameters theta 
                             Pr_Y0_theta <- function(my_data){
                               my_Pr_Y0_theta <- (1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data) + (1-Pr_Y1_Lambda0(my_data$epsilon))*(1-Pr_Lambda1(my_data))
                               return(my_Pr_Y0_theta)
                             }
                             
                             ### Function for calculating p_i 
                             calculate_p <- function(my_data){
                               my_data %>% 
                                 dplyr::mutate(p = ifelse(Y == 1,(Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data))/ Pr_Y1_theta(my_data), 
                                                          ((1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data))/Pr_Y0_theta(my_data))) -> my_p_table
                               return(my_p_table)
                             }
                             
                             ### Function for calculating the complete log-likelihood  
                             e_step_LL <- function(my_data){
                               Y1 <- my_data %>% dplyr::filter(Y == 1)
                               Y0 <- my_data %>% dplyr::filter(Y == 0)
                               my_data %>% 
                                 dplyr::select(beta0,beta1) %>% 
                                 dplyr::distinct() %>% 
                                 as.numeric() -> my_beta
                               Y1_X1 <- cbind(1, as.matrix(Y1['X1']))
                               Y0_X1 <- cbind(1, as.matrix(Y0['X1']))
                               
                               LL_Y1 <- sum((Y1$p*(Y1$Y*Y1$fMi - log(1 + exp(Y1$fMi)))) + 
                                              ((1-Y1$p)*(Y1$Y*Y1$epsilon - log(1+exp(Y1$epsilon)))) + 
                                              (Y1$p*(Y1_X1 %*% my_beta) - log(1+exp(Y1_X1 %*% my_beta)))) 
                               
                               LL_Y0 <- sum((Y0$p)*(Y0$Y*Y0$fMi - log(1 + exp(Y0$fMi))) + 
                                              ((1-Y0$p)*(Y0$Y*Y0$epsilon - log(1+exp(Y0$epsilon)))) + 
                                              (Y0$p*Y0_X1 %*% my_beta - log(1+exp(Y0_X1 %*% my_beta)))) 
                               
                               Expected_LL = LL_Y1 + LL_Y0   
                               return(Expected_LL)
                             }
                             
                             ### Function for Maximization of vBeta
                             vBeta <- function(my_data){
                               ### logistic regression:
                               # my_updated_beta <- glm(p~X1, data = my_data, family = binomial)
                               ### Firth-penalized logistic regression:
                               my_updated_beta <- logistf::logistf(p~X1, data = my_data)
                               my_vBeta <-c(my_updated_beta$coefficients[[1]],my_updated_beta$coefficients[[2]])
                               return(my_vBeta)
                             }
                             
                             ### Function for Maximization of fMi 
                             vfMi <- function(my_data){
                               my_data %>% dplyr::arrange(coverage) -> my_data
                               my_data$Y -> y 
                               my_data$p -> w1
                               n <- nrow(my_data)
                               
                               Atot <- cbind(1:(n-1),2:n) # define monotonicity 
                               
                               fit <- isotone::activeSet(isomat = Atot,
                                                         mySolver = fSolver,
                                                         fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/50)^2)),
                                                         gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (x/25)*sinh((x/50)^2),
                                                         y = y,
                                                         weights = w1)
                               # a =  25
                               # fit <- isotone::activeSet(isomat = Atot,
                               #                           mySolver = fSolver,
                               #                           fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/a)^4)),
                               #                           gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (4*(x^3)/a)*sinh((x/50)^4),
                               #                           y = y,
                               #                           weights = w1)
                               fit2 <- cbind(fit$y,fit$x,my_data$coverage) %>% as.data.frame()
                               colnames(fit2) <- c("Y","fMi","coverage")
                               return(fit2)
                             }
                             
                             num <- nrow(draw)
                             draw$fMi <- rep(1,num) # I'm going to generate some random fMi's to initialize
                             draw$epsilon <- 0
                             # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                             
                             #t <- 100 # how many iterations for algorithm to converge
                             my_data1 <- calculate_p(draw) # initial dataset calculation of p's 
                             my_data1 <- my_data1 %>% 
                               dplyr::mutate(iteration = 1)
                             LL_list <- list() # List to store LL's 
                             
                             mydata_list <- list() # List to store datasets corresponding to LL's 
                             
                             i <- 1 
                             rel_tol <- 1
                             #abs_tol <- 1
                             LL_list[i] <- e_step_LL(eval(as.symbol(paste0("my_data",i))))
                             mydata_list[i] <- list(eval(as.symbol(paste0("my_data",i))))
                             
                             while (rel_tol > 0.00001) {
                               j <- i+1 # increment by i + 1 
                               # I need to calculate the p's which I did outside this loop for the initial dataset
                               # Then calculate my LL for that iteration 
                               
                               # Using p's I can update my parameters for the next dataset j to use 
                               assign(paste0("nextbetas_",j), vBeta(eval(as.symbol(paste0("my_data",i)))))
                               assign(paste0("nextfMis_",j), vfMi(eval(as.symbol(paste0("my_data",i)))))
                               
                               # Then I need to Drop the previous variables that will be replaced with new ones
                               drop.cols <- c('fMi','beta0','beta1','p')
                               assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",i))) %>% dplyr::select(-drop.cols))
                               # Then Add the new betas and fMis to the dataframe 
                               assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",j))) %>% 
                                        dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetas_",j)))[1], 
                                                      beta1 = eval(as.symbol(paste0("nextbetas_",j)))[2], 
                                                      iteration = j))
                               mergeCols <- c("coverage","Y")
                               assign(paste0("my_data",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMis_",j))), eval(as.symbol(paste0("my_data",j))), by = mergeCols))
                               #Now I have an updated dataframe of fMis and betas 
                               
                               #Next I need to Calculate the new P's for this 
                               
                               assign(paste0("my_data",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_data",j)))))
                               
                               mydata_list[j] <- list(eval(as.symbol(paste0("my_data",j))))
                               # Then calculate my LL for that iteration 
                               LL_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                               
                               abs_tol <- abs(LL_list[[j]] - LL_list[[i]])
                               rel_tol <- (abs(LL_list[[j]] - LL_list[[i]]))/abs(LL_list[[i]])
                               
                               i = i + 1
                               print(j)
                               if (j == 1000){
                                 rel_tol <- -1 
                                 LL_list[[j]] <- "max iterations reached"
                               }
                             }
                             length(mydata_list) -> em_iteration
                             
                             mydata_list[[em_iteration]] %>% 
                               dplyr::select(beta0,beta1) %>% 
                               dplyr::distinct() %>% 
                               as.numeric() -> my_em_betas
                             
                             
                             output <- list("beta0" = my_em_betas[[1]],
                                            "beta1" = my_em_betas[[2]], 
                                            "LL" = LL_list[[em_iteration]], 
                                            "iterations" = em_iteration)
                           }
)

## old method
david_method_nopenalty <- new_method("david_method_nopenalty","david_method_nopenalty",
                                     method = function(model,draw){
                                       ### Functions for the 3 hierarchical model components 
                                       # Probability of Y = 1 given Lambda = 1
                                       library(isotone)
                                       Pr_Y1_Lambda1 <- function(my_fMi) {
                                         my_Pr_Y1_Lambda1 <- exp(my_fMi)/(1+exp(my_fMi))
                                         return(my_Pr_Y1_Lambda1)
                                       }
                                       # Probability of Lambda = 1 
                                       Pr_Lambda1 <- function(my_data) {
                                         mX1 = as.matrix(my_data['X1'])
                                         mX1 = cbind(1, mX1)
                                         my_data %>% 
                                           dplyr::select(beta0,beta1) %>% 
                                           dplyr::distinct() %>% 
                                           as.numeric() -> my_beta
                                         my_Pr_Lambda1 <- exp(mX1 %*% my_beta)/(1+exp(mX1 %*% my_beta))
                                         return(my_Pr_Lambda1)
                                       }
                                       # Probability of Y = 1 given Lambda = 0 
                                       Pr_Y1_Lambda0 <- function(my_epsilon) {
                                         #my_Pr_Y1_Lambda0 <- exp(my_epsilon)/(1+exp(my_epsilon))
                                         # for now I'm going to set this parameterization to 0 
                                         my_Pr_Y1_Lambda0 <- 0 
                                         return(my_Pr_Y1_Lambda0)
                                       }
                                       
                                       ### Functions for calculating the denominators of p_i 
                                       # Prob of Y = 1 given parameters theta 
                                       Pr_Y1_theta <- function(my_data){
                                         my_Pr_Y1_theta <- Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data) + Pr_Y1_Lambda0(my_data$epsilon)*(1-Pr_Lambda1(my_data))
                                         return(my_Pr_Y1_theta)
                                       }
                                       #Prob of Y = 0 given parameters theta 
                                       Pr_Y0_theta <- function(my_data){
                                         my_Pr_Y0_theta <- (1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data) + (1-Pr_Y1_Lambda0(my_data$epsilon))*(1-Pr_Lambda1(my_data))
                                         return(my_Pr_Y0_theta)
                                       }
                                       
                                       ### Function for calculating p_i 
                                       calculate_p <- function(my_data){
                                         my_data %>% 
                                           dplyr::mutate(p = ifelse(Y == 1,(Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data))/ Pr_Y1_theta(my_data), 
                                                                    ((1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data))/Pr_Y0_theta(my_data))) -> my_p_table
                                         return(my_p_table)
                                       }
                                       
                                       ### Function for calculating the complete log-likelihood  
                                       e_step_LL <- function(my_data){
                                         Y1 <- my_data %>% dplyr::filter(Y == 1)
                                         Y0 <- my_data %>% dplyr::filter(Y == 0)
                                         my_data %>% 
                                           dplyr::select(beta0,beta1) %>% 
                                           dplyr::distinct() %>% 
                                           as.numeric() -> my_beta
                                         Y1_X1 <- cbind(1, as.matrix(Y1['X1']))
                                         Y0_X1 <- cbind(1, as.matrix(Y0['X1']))
                                         
                                         LL_Y1 <- sum((Y1$p*(Y1$Y*Y1$fMi - log(1 + exp(Y1$fMi)))) + 
                                                        ((1-Y1$p)*(Y1$Y*Y1$epsilon - log(1+exp(Y1$epsilon)))) + 
                                                        (Y1$p*(Y1_X1 %*% my_beta) - log(1+exp(Y1_X1 %*% my_beta)))) 
                                         
                                         LL_Y0 <- sum((Y0$p)*(Y0$Y*Y0$fMi - log(1 + exp(Y0$fMi))) + 
                                                        ((1-Y0$p)*(Y0$Y*Y0$epsilon - log(1+exp(Y0$epsilon)))) + 
                                                        (Y0$p*Y0_X1 %*% my_beta - log(1+exp(Y0_X1 %*% my_beta)))) 
                                         
                                         Expected_LL = LL_Y1 + LL_Y0   
                                         return(Expected_LL)
                                       }
                                       
                                       ### Function for Maximization of vBeta
                                       vBeta <- function(my_data){
                                         ### logistic regression:
                                         # my_updated_beta <- glm(p~X1, data = my_data, family = binomial)
                                         ### Firth-penalized logistic regression:
                                         my_updated_beta <- logistf::logistf(p~X1, data = my_data)
                                         my_vBeta <-c(my_updated_beta$coefficients[[1]],my_updated_beta$coefficients[[2]])
                                         return(my_vBeta)
                                       }
                                       
                                       ### Function for Maximization of fMi 
                                       vfMi <- function(my_data){
                                         my_data %>% dplyr::arrange(coverage) -> my_data
                                         my_data$Y -> y 
                                         my_data$p -> w1
                                         n <- nrow(my_data)
                                         
                                         Atot <- cbind(1:(n-1),2:n) # define monotonicity 
                                         fit <- isotone::activeSet(isomat = Atot, mySolver = fSolver, fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))), 
                                                                   gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))), y = y, weights = w1)
                                         # a =  25
                                         # fit <- isotone::activeSet(isomat = Atot,
                                         #                           mySolver = fSolver,
                                         #                           fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/a)^4)),
                                         #                           gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (4*(x^3)/a)*sinh((x/50)^4),
                                         #                           y = y,
                                         #                           weights = w1)
                                         fit2 <- cbind(fit$y,fit$x,my_data$coverage) %>% as.data.frame()
                                         colnames(fit2) <- c("Y","fMi","coverage")
                                         return(fit2)
                                       }
                                       
                                       num <- nrow(draw)
                                       draw$fMi <- rep(1,num) # I'm going to generate some random fMi's to initialize
                                       draw$epsilon <- 0
                                       # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                       
                                       #t <- 100 # how many iterations for algorithm to converge
                                       my_data1 <- calculate_p(draw) # initial dataset calculation of p's 
                                       my_data1 <- my_data1 %>% 
                                         dplyr::mutate(iteration = 1)
                                       LL_list <- list() # List to store LL's 
                                       
                                       mydata_list <- list() # List to store datasets corresponding to LL's 
                                       
                                       i <- 1 
                                       rel_tol <- 1
                                       #abs_tol <- 1
                                       LL_list[i] <- e_step_LL(eval(as.symbol(paste0("my_data",i))))
                                       mydata_list[i] <- list(eval(as.symbol(paste0("my_data",i))))
                                       
                                       while (rel_tol > 0.00001) {
                                         j <- i+1 # increment by i + 1 
                                         # I need to calculate the p's which I did outside this loop for the initial dataset
                                         # Then calculate my LL for that iteration 
                                         
                                         # Using p's I can update my parameters for the next dataset j to use 
                                         assign(paste0("nextbetas_",j), vBeta(eval(as.symbol(paste0("my_data",i)))))
                                         assign(paste0("nextfMis_",j), vfMi(eval(as.symbol(paste0("my_data",i)))))
                                         
                                         # Then I need to Drop the previous variables that will be replaced with new ones
                                         drop.cols <- c('fMi','beta0','beta1','p')
                                         assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",i))) %>% dplyr::select(-drop.cols))
                                         # Then Add the new betas and fMis to the dataframe 
                                         assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",j))) %>% 
                                                  dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetas_",j)))[1], 
                                                                beta1 = eval(as.symbol(paste0("nextbetas_",j)))[2], 
                                                                iteration = j))
                                         mergeCols <- c("coverage","Y")
                                         assign(paste0("my_data",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMis_",j))), eval(as.symbol(paste0("my_data",j))), by = mergeCols))
                                         #Now I have an updated dataframe of fMis and betas 
                                         
                                         #Next I need to Calculate the new P's for this 
                                         
                                         assign(paste0("my_data",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_data",j)))))
                                         
                                         mydata_list[j] <- list(eval(as.symbol(paste0("my_data",j))))
                                         # Then calculate my LL for that iteration 
                                         LL_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                         
                                         abs_tol <- abs(LL_list[[j]] - LL_list[[i]])
                                         rel_tol <- (abs(LL_list[[j]] - LL_list[[i]]))/abs(LL_list[[i]])
                                         
                                         i = i + 1
                                         print(j)
                                         if (j == 1000){
                                           rel_tol <- -1 
                                           LL_list[[j]] <- "max iterations reached"
                                         }
                                       }
                                       length(mydata_list) -> em_iteration
                                       
                                       mydata_list[[em_iteration]] %>% 
                                         dplyr::select(beta0,beta1) %>% 
                                         dplyr::distinct() %>% 
                                         as.numeric() -> my_em_betas
                                       
                                       
                                       output <- list("beta0" = my_em_betas[[1]],
                                                      "beta1" = my_em_betas[[2]], 
                                                      "LL" = LL_list[[em_iteration]], 
                                                      "iterations" = em_iteration)
                                     }
)


david_method2 <- new_method("david_method2","david_method2",
                            method = function(model,draw){
                              ### Functions for the 3 hierarchical model components 
                              # Probability of Y = 1 given Lambda = 1
                              library(isotone)
                              
                              Pr_Y1_Lambda1 <- function(my_fMi) {
                                my_Pr_Y1_Lambda1 <- exp(my_fMi)/(1+exp(my_fMi))
                                return(my_Pr_Y1_Lambda1)
                              }
                              # Probability of Lambda = 1 
                              
                              Pr_Lambda1 <- function(my_data) {
                                mX1 = as.matrix(my_data['X1'])
                                mX1 = cbind(1, mX1)
                                my_data %>% 
                                  dplyr::select(beta0,beta1) %>% 
                                  dplyr::distinct() %>% 
                                  as.numeric() -> my_beta
                                my_Pr_Lambda1 <- exp(mX1 %*% my_beta)/(1+exp(mX1 %*% my_beta))
                                return(my_Pr_Lambda1)
                              }
                              
                              Pr_Lambda1_null <- function(my_data) {
                                my_data %>% 
                                  dplyr::select(beta0) %>% 
                                  dplyr::distinct() %>% 
                                  as.numeric() -> my_beta
                                mX1 = matrix(1,nrow(my_data),1)
                                my_Pr_Lambda1_null <- exp(mX1 %*% my_beta)/(1+exp(mX1 %*% my_beta))
                                return(my_Pr_Lambda1_null)
                              }
                              
                              # Probability of Y = 1 given Lambda = 0 
                              Pr_Y1_Lambda0 <- function(my_epsilon) {
                                #my_Pr_Y1_Lambda0 <- exp(my_epsilon)/(1+exp(my_epsilon))
                                # for now I'm going to set this parameterization to 0 
                                my_Pr_Y1_Lambda0 <- 0 
                                return(my_Pr_Y1_Lambda0)
                              }
                              
                              ### Functions for calculating the denominators of p_i 
                              # Prob of Y = 1 given parameters theta (DENOMINATOR of p_i)
                              Pr_Y1_theta <- function(my_data){
                                my_Pr_Y1_theta <- Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data) + Pr_Y1_Lambda0(my_data$epsilon)*(1-Pr_Lambda1(my_data))
                                return(my_Pr_Y1_theta)
                              }
                              #Prob of Y = 0 given parameters theta (DENOMINATOR of p_i)
                              Pr_Y0_theta <- function(my_data){
                                my_Pr_Y0_theta <- (1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data) + (1-Pr_Y1_Lambda0(my_data$epsilon))*(1-Pr_Lambda1(my_data))
                                return(my_Pr_Y0_theta)
                              }
                              
                              ### Function for calculating p_i 
                              calculate_p <- function(my_data){
                                my_data %>% 
                                  dplyr::mutate(p = ifelse(Y == 1,(Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data))/ Pr_Y1_theta(my_data), 
                                                           ((1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data))/Pr_Y0_theta(my_data))) -> my_p_table
                                return(my_p_table)
                              }
                              
                              ### Function for calculating the expected complete log-likelihood  
                              e_step_LL <- function(my_data){
                                Y1 <- my_data %>% dplyr::filter(Y == 1)
                                Y0 <- my_data %>% dplyr::filter(Y == 0)
                                my_data %>% 
                                  dplyr::select(beta0,beta1) %>% 
                                  dplyr::distinct() %>% 
                                  as.numeric() -> my_beta
                                Y1_X1 <- cbind(1, as.matrix(Y1['X1']))
                                Y0_X1 <- cbind(1, as.matrix(Y0['X1']))
                                
                                LL_Y1 <- sum((Y1$p*(Y1$Y*Y1$fMi - log(1 + exp(Y1$fMi)))) + 
                                               ((1-Y1$p)*(Y1$Y*Y1$epsilon - log(1+exp(Y1$epsilon)))) + 
                                               (Y1$p*(Y1_X1 %*% my_beta) - log(1+exp(Y1_X1 %*% my_beta)))) 
                                
                                LL_Y0 <- sum((Y0$p)*(Y0$Y*Y0$fMi - log(1 + exp(Y0$fMi))) + 
                                               ((1-Y0$p)*(Y0$Y*Y0$epsilon - log(1+exp(Y0$epsilon)))) + 
                                               (Y0$p*Y0_X1 %*% my_beta - log(1+exp(Y0_X1 %*% my_beta)))) 
                                
                                Expected_LL = LL_Y1 + LL_Y0   
                                return(Expected_LL)
                              }
                              
                              ### Function for calculating the incomplete log-likelihood 
                              incomplete_LL <- function(par, my_data){
                                e = 0 # e is parameter argument 1 
                                vBeta = c(par[1],par[2]) # my beta vector is parameter arguments 2 and 3 
                                Y0 <- my_data %>% dplyr::filter(Y == 0)
                                Y1 <- my_data %>% dplyr::filter(Y == 1)
                                
                                Mi0 = as.matrix(Y0['fMi']) # f_tilde for Y=0
                                # Mi0 = log(Mi0)
                                # add an intercept to the predictor variables
                                #Mi0 = cbind(1, Mi0)
                                mX0 = as.matrix(Y0['X1'])
                                mX0 = cbind(1, mX0)
                                
                                Mi1 = as.matrix(Y1['fMi']) # f_tilde for Y = 1
                                # Mi1 = log(Mi1) # modeling the log coverage
                                #Mi1 = cbind(1, Mi1)
                                mX1 = as.matrix(Y1['X1'])
                                mX1 = cbind(1, mX1)
                                
                                LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                (1-(exp(Mi0)/(1+exp(Mi0))))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                  sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                            (1/(1+exp(-(mX1 %*% vBeta))))*(exp(Mi1)/(1+exp(Mi1)))))
                                return(LL)
                              }
                              
                              ### Function for Maximization of vBeta
                              vBeta <- function(my_data){
                                ### logistic regression:
                                # my_updated_beta <- glm(p~X1, data = my_data, family = binomial)
                                ### Firth-penalized logistic regression:
                                my_updated_beta <- logistf::logistf(p~X1, data = my_data)
                                my_vBeta <-c(my_updated_beta$coefficients[[1]],my_updated_beta$coefficients[[2]])
                                return(my_vBeta)
                              }
                              
                              ### Function for maximization of vBeta Null model 
                              vBeta_null <- function(my_data){
                                ### logistic regression:
                                # my_updated_beta <- glm(p~X1, data = my_data, family = binomial)
                                ### Firth-penalized logistic regression:
                                my_updated_beta <- logistf::logistf(p~1, data = my_data)
                                my_vBeta <- my_updated_beta$coefficients[[1]]
                                return(my_vBeta)
                              }
                              
                              ### Function for Maximization of fMi 
                              vfMi <- function(my_data){
                                my_data %>% dplyr::arrange(coverage) -> my_data
                                my_data$Y -> y 
                                my_data$p -> w1
                                n <- nrow(my_data)
                                
                                Atot <- cbind(1:(n-1),2:n) # define monotonicity 
                                
                                fit <- isotone::activeSet(isomat = Atot,
                                                          mySolver = fSolver,
                                                          fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/50)^2)),
                                                          gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (x/25)*sinh((x/50)^2),
                                                          y = y,
                                                          weights = w1)
                                # a =  25
                                # fit <- isotone::activeSet(isomat = Atot,
                                #                           mySolver = fSolver,
                                #                           fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/a)^4)),
                                #                           gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (4*(x^3)/a)*sinh((x/50)^4),
                                #                           y = y,
                                #                           weights = w1)
                                fit2 <- cbind(fit$y,fit$x,my_data$coverage) %>% as.data.frame()
                                colnames(fit2) <- c("Y","fMi","coverage")
                                return(fit2)
                              }
                              
                              if(length(draw$Y %>% unique()) == 1){
                                my_results <- c(0,0,0,0,0)
                                my_pvalue <- 1
                                my_LRT <- 0
                              } else if(length(draw$X1 %>% unique()) == 1) {
                                my_results <- c(0,0,0,0,0)
                                my_pvalue <- 1
                                my_LRT <- 0
                                
                              } else {
                                
                                num <- nrow(draw)
                                draw$epsilon <- 0
                                
                                init_list <- list()
                                init_list[[1]] <- rep(1,num)
                                #init_list[[2]] <- seq(-1,1,length = num)
                                
                                final_LL_list <- list() 
                                
                                for (k in 1:length(init_list)){  
                                  
                                  draw %>% dplyr::arrange(coverage) -> draw
                                  draw$fMi <- init_list[[k]] # Start with a uniform fMi's
                                  # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                  
                                  #t <- 100 # how many iterations for algorithm to converge
                                  my_data1 <- calculate_p(draw) # initial dataset calculation of p's 
                                  my_data1 <- my_data1 %>% 
                                    dplyr::mutate(iteration = 1)
                                  
                                  LL_list <- list() # List to store LL's 
                                  
                                  mydata_list <- list() # List to store datasets corresponding to LL's 
                                  
                                  ftilde_list <- list() 
                                  
                                  i <- 1 
                                  beta1_change <- 1
                                  beta0_change <- 1
                                  fmi_change <- 1
                                  rel_tol <- 1
                                  #abs_tol <- 1
                                  eval(as.symbol(paste0("my_data",i))) %>% 
                                    dplyr::select(beta0,beta1) %>% 
                                    dplyr::distinct() %>% 
                                    as.numeric() -> my_init_betas
                                  
                                  LL_list[i] <- incomplete_LL(my_init_betas,eval(as.symbol(paste0("my_data",i))))
                                  mydata_list[i] <- list(eval(as.symbol(paste0("my_data",i))))
                                  ftilde_list[i] <- -1 
                                  while (beta1_change > 0.0001 | beta0_change > 0.0001 | fmi_change > 0.01) {
                                    j <- i+1 # increment by i + 1 
                                    # I need to calculate the p's which I did outside this loop for the initial dataset
                                    # Then calculate my LL for that iteration 
                                    
                                    # Using p's I can update my parameters for the next dataset j to use 
                                    assign(paste0("nextbetas_",j), vBeta(eval(as.symbol(paste0("my_data",i)))))
                                    assign(paste0("nextfMis_",j), vfMi(eval(as.symbol(paste0("my_data",i)))))
                                    
                                    # Then I need to Drop the previous variables that will be replaced with new ones
                                    drop.cols <- c('fMi','beta0','beta1','p')
                                    assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",i))) %>% dplyr::select(-drop.cols))
                                    # Then Add the new betas and fMis to the dataframe 
                                    assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",j))) %>% 
                                             dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetas_",j)))[1], 
                                                           beta1 = eval(as.symbol(paste0("nextbetas_",j)))[2], 
                                                           iteration = j))
                                    mergeCols <- c("coverage","Y")
                                    assign(paste0("my_data",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMis_",j))), eval(as.symbol(paste0("my_data",j))), by = mergeCols))
                                    #Now I have an updated dataframe of fMis and betas 
                                    
                                    #Next I need to Calculate the new P's for this 
                                    
                                    assign(paste0("my_data",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_data",j)))))
                                    
                                    mydata_list[j] <- list(eval(as.symbol(paste0("my_data",j))))
                                    # Then calculate my LL for that iteration 
                                    eval(as.symbol(paste0("my_data",j))) %>% 
                                      dplyr::select(beta0,beta1) %>% 
                                      dplyr::distinct() %>% 
                                      as.numeric() -> my_j_betas
                                    
                                    LL_list[j] <- incomplete_LL(my_j_betas,eval(as.symbol(paste0("my_data",j))))
                                    
                                    #  LL_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                    
                                    abs_tol <- abs(LL_list[[j]] - LL_list[[i]])
                                    rel_tol <- (abs(LL_list[[j]] - LL_list[[i]]))/abs(LL_list[[i]])
                                    
                                    betas_next <- mydata_list[[j]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                    betas_previous <- mydata_list[[i]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                    
                                    beta0_change <- abs(betas_next[[1]]-betas_previous[[1]])
                                    beta1_change <- abs(betas_next[[2]]-betas_previous[[2]])
                                    
                                    fmis_next <- mydata_list[[j]] %>% dplyr::select(fMi) 
                                    fmis_previous <- mydata_list[[i]] %>% dplyr::select(fMi)
                                    
                                    fmis_intermediate <- abs(fmis_next-fmis_previous)/abs(fmis_previous)
                                    max(fmis_intermediate)
                                    fmi_change <- max(fmis_intermediate)
                                    
                                    ftilde_list[j] <- fmi_change
                                    
                                    i = i + 1
                                    
                                    print(j)
                                    if (j == 500){
                                      rel_tol <- -1 
                                      beta1_change <- -1
                                      beta0_change <- -1
                                      fmi_change <- -1 
                                      #LL_list[[j]] <- "max iterations reached"
                                      # do.call("rbind", ftilde_list) -> ftilde_table
                                      # ftilde_table %>% 
                                      #   dplyr::as_data_frame() %>% 
                                      #   dplyr::mutate(iteration = 1:dplyr::n()) %>% 
                                      #   dplyr::filter(iteration > 1) -> ftilde_change_plot
                                      # ftilde_change_plot[which.min(ftilde_change_plot$V1),][[1]] -> ftilde_list[[j]]
                                    }
                                  }
                                  length(mydata_list) -> em_iteration
                                  
                                  mydata_list[[em_iteration]] %>% 
                                    dplyr::select(beta0,beta1) %>% 
                                    dplyr::distinct() %>% 
                                    as.numeric() -> my_em_betas
                                  
                                  cbind(my_em_betas[[1]], my_em_betas[[2]],LL_list[[em_iteration]], 
                                        em_iteration[[1]],ftilde_list[[em_iteration]]) -> final_LL_list[[k]]
                                }
                                
                                do.call(rbind, final_LL_list) %>% 
                                  dplyr::as_tibble() %>% 
                                  dplyr::filter(!is.na(V1)) %>% 
                                  dplyr::distinct() -> results_tibble # but really there's only one 
                                results_tibble[which.max(results_tibble$V5),] -> my_results 
                                if (dim(my_results)[[1]] == 0){
                                  print("unrestricted model didn't converge")
                                  my_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                } else {
                                  print("unrestricted model completed!")
                                }
                                
                                ############# for the restricted/null model ###############
                                
                                num <- nrow(draw)
                                draw$epsilon <- 0
                                
                                init_list <- list()
                                init_list[[1]] <- rep(1,num)
                                #init_list[[2]] <- seq(-1,1,length = num) 
                                ## previously when we had multiple initial starts
                                
                                restricted_LL_list <- list() 
                                
                                for (k in 1:length(init_list)){  
                                  
                                  draw %>% dplyr::arrange(coverage) -> draw
                                  draw$fMi <- init_list[[k]] # Start with a uniform fMi's
                                  # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                  draw %>% 
                                    dplyr::select(-beta1) %>% 
                                    dplyr::mutate(beta1 = 0) -> null_draw
                                  #t <- 100 # how many iterations for algorithm to converge
                                  my_nulldata1 <- calculate_p(null_draw) # initial dataset calculation of p's 
                                  my_nulldata1 <- my_nulldata1 %>% 
                                    dplyr::mutate(iteration = 1)
                                  
                                  LL_restricted_list <- list() # List to store LL's 
                                  
                                  mydata_restricted_list <- list() # List to store datasets corresponding to LL's 
                                  
                                  ftilde_restricted_list <- list() 
                                  
                                  i <- 1 
                                  beta0_null_change <- 1
                                  fmi_null_change <- 1
                                  #rel_null_tol <- 1
                                  #abs_tol <- 1
                                  eval(as.symbol(paste0("my_nulldata",i))) %>% 
                                    dplyr::select(beta0,beta1) %>% 
                                    dplyr::distinct() %>% 
                                    as.numeric()  -> my_init_betas
                                  
                                  LL_restricted_list[i] <- incomplete_LL(my_init_betas,eval(as.symbol(paste0("my_nulldata",i))))
                                  
                                  mydata_restricted_list[i] <- list(eval(as.symbol(paste0("my_nulldata",i))))
                                  ftilde_restricted_list[i] <- -1 
                                  while (beta0_null_change > 0.0001 | fmi_null_change > 0.01) {
                                    j <- i+1 # increment by i + 1 
                                    # I need to calculate the p's which I did outside this loop for the initial dataset
                                    # Then calculate my LL for that iteration 
                                    
                                    # Using p's I can update my parameters for the next dataset j to use 
                                    assign(paste0("nextbetasnull_",j), vBeta_null(eval(as.symbol(paste0("my_nulldata",i)))))
                                    assign(paste0("nextfMisnull_",j), vfMi(eval(as.symbol(paste0("my_nulldata",i)))))
                                    
                                    # Then I need to Drop the previous variables that will be replaced with new ones
                                    drop.cols <- c('fMi','beta0','beta1','p')
                                    assign(paste0("my_nulldata",j), eval(as.symbol(paste0("my_nulldata",i))) %>% dplyr::select(-drop.cols))
                                    # Then Add the new betas and fMis to the dataframe 
                                    assign(paste0("my_nulldata",j), eval(as.symbol(paste0("my_nulldata",j))) %>% 
                                             dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetasnull_",j)))[1],
                                                           beta1 = 0 , 
                                                           iteration = j))
                                    mergeCols <- c("coverage","Y")
                                    assign(paste0("my_nulldata",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMisnull_",j))), eval(as.symbol(paste0("my_nulldata",j))), by = mergeCols))
                                    #Now I have an updated dataframe of fMis and betas 
                                    
                                    #Next I need to Calculate the new P's for this 
                                    
                                    assign(paste0("my_nulldata",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_nulldata",j)))))
                                    
                                    mydata_restricted_list[j] <- list(eval(as.symbol(paste0("my_nulldata",j))))
                                    # Then calculate my LL for that iteration 
                                    eval(as.symbol(paste0("my_nulldata",j))) %>% 
                                      dplyr::select(beta0,beta1) %>% 
                                      dplyr::distinct() %>% 
                                      as.numeric() -> my_j_betas 
                                    
                                    LL_restricted_list[j] <- incomplete_LL(my_j_betas,eval(as.symbol(paste0("my_nulldata",i))))
                                    
                                    #LL_restricted_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                    
                                    abs_tol <- abs(LL_restricted_list[[j]] - LL_restricted_list[[i]])
                                    rel_tol <- (abs(LL_restricted_list[[j]] - LL_restricted_list[[i]]))/abs(LL_restricted_list[[i]])
                                    
                                    betas_null_next <- mydata_restricted_list[[j]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                    betas_null_previous <- mydata_restricted_list[[i]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                    
                                    beta0_null_change <- abs(betas_null_next[[1]]-betas_null_previous[[1]])
                                    beta1_null_change <- abs(betas_null_next[[2]]-betas_null_previous[[2]])
                                    
                                    fmis_null_next <- mydata_restricted_list[[j]] %>% dplyr::select(fMi) 
                                    fmis_null_previous <- mydata_restricted_list[[i]] %>% dplyr::select(fMi)
                                    
                                    fmis_null_intermediate <- abs(fmis_null_next-fmis_null_previous)/abs(fmis_null_previous)
                                    max(fmis_null_intermediate)
                                    fmi_null_change <- max(fmis_null_intermediate)
                                    
                                    ftilde_restricted_list[j] <- fmi_null_change
                                    
                                    i = i + 1
                                    
                                    print(j)
                                    if (j == 500){
                                      rel_tol <- -1 
                                      beta1_null_change <- -1
                                      beta0_null_change <- -1
                                      fmi_null_change <- -1 
                                      #LL_list[[j]] <- "max iterations reached"
                                      # do.call("rbind", ftilde_restricted_list) -> ftilde_table
                                      # ftilde_table %>% 
                                      #   dplyr::as_data_frame() %>% 
                                      #   dplyr::mutate(iteration = 1:dplyr::n()) %>% 
                                      #   dplyr::filter(iteration > 1) -> ftilde_change_plot
                                      # ftilde_change_plot[which.min(ftilde_change_plot$V1),][[1]] -> ftilde_list[[j]]
                                    }
                                  }
                                  
                                  length(mydata_restricted_list) -> em_restricted_iteration
                                  
                                  mydata_restricted_list[[em_restricted_iteration]] %>% 
                                    dplyr::select(beta0,beta1) %>% 
                                    dplyr::distinct() %>% 
                                    as.numeric() -> my_em_restricted_betas
                                  
                                  cbind(my_em_restricted_betas[[1]], my_em_restricted_betas[[2]],LL_restricted_list[[em_restricted_iteration]], 
                                        em_restricted_iteration[[1]],ftilde_restricted_list[[em_restricted_iteration]]) -> restricted_LL_list[[k]]
                                }
                                
                                do.call(rbind, restricted_LL_list) %>% 
                                  dplyr::as_tibble() %>% 
                                  dplyr::filter(!is.na(V1)) %>% 
                                  dplyr::distinct() -> restricted_results_tibble
                                
                                restricted_results_tibble[which.max(restricted_results_tibble$V5),] -> my_restricted_results 
                                ## ^^ from when we had multiple initial starts and were trying to use the highest LL but in this case we only have 1 start
                                
                                if (dim(my_restricted_results)[[1]] == 0){
                                  print("restricted model didn't converge")
                                  my_restricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                } else {
                                  print("restricted model completed!")
                                }
                                
                                my_LRT = 2*abs(my_results$V3 - my_restricted_results$V3)
                                
                                my_pvalue <- stats::pchisq(my_LRT, df = 1, lower.tail = FALSE)
                                
                              }
                              
                              output <- list("beta0" = my_results[[1]],
                                             "beta1" = my_results[[2]], 
                                             "LRT" = my_LRT, 
                                             "iterations" = my_results[[4]], 
                                             "ftilde_change" = my_results[[5]], 
                                             "pvalue"= my_pvalue)
                            })


em_method_wasserman <- new_method("em_method_wasserman","em_method_wasserman",
                                  method = function(model,draw){
                                    ### Functions for the 3 hierarchical model components 
                                    # Probability of Y = 1 given Lambda = 1
                                    library(isotone)
                                    Pr_Y1_Lambda1 <- function(my_fMi) {
                                      my_Pr_Y1_Lambda1 <- exp(my_fMi)/(1+exp(my_fMi))
                                      return(my_Pr_Y1_Lambda1)
                                    }
                                    # Probability of Lambda = 1 
                                    Pr_Lambda1 <- function(my_data) {
                                      mX1 = as.matrix(my_data['X1'])
                                      mX1 = cbind(1, mX1)
                                      my_data %>% 
                                        dplyr::select(beta0,beta1) %>% 
                                        dplyr::distinct() %>% 
                                        as.numeric() -> my_beta
                                      my_Pr_Lambda1 <- exp(mX1 %*% my_beta)/(1+exp(mX1 %*% my_beta))
                                      return(my_Pr_Lambda1)
                                    }
                                    
                                    Pr_Lambda1_null <- function(my_data) {
                                      my_data %>% 
                                        dplyr::select(beta0) %>% 
                                        dplyr::distinct() %>% 
                                        as.numeric() -> my_beta
                                      mX1 = matrix(1,nrow(my_data),1)
                                      my_Pr_Lambda1_null <- exp(mX1 %*% my_beta)/(1+exp(mX1 %*% my_beta))
                                      return(my_Pr_Lambda1_null)
                                    }
                                    
                                    # Probability of Y = 1 given Lambda = 0 
                                    Pr_Y1_Lambda0 <- function(my_epsilon) {
                                      #my_Pr_Y1_Lambda0 <- exp(my_epsilon)/(1+exp(my_epsilon))
                                      # for now I'm going to set this parameterization to 0 
                                      my_Pr_Y1_Lambda0 <- 0 
                                      return(my_Pr_Y1_Lambda0)
                                    }
                                    
                                    ### Functions for calculating the denominators of p_i 
                                    # Prob of Y = 1 given parameters theta (DENOMINATOR of p_i)
                                    Pr_Y1_theta <- function(my_data){
                                      my_Pr_Y1_theta <- Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data) + Pr_Y1_Lambda0(my_data$epsilon)*(1-Pr_Lambda1(my_data))
                                      return(my_Pr_Y1_theta)
                                    }
                                    #Prob of Y = 0 given parameters theta (DENOMINATOR of p_i)
                                    Pr_Y0_theta <- function(my_data){
                                      my_Pr_Y0_theta <- (1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data) + (1-Pr_Y1_Lambda0(my_data$epsilon))*(1-Pr_Lambda1(my_data))
                                      return(my_Pr_Y0_theta)
                                    }
                                    
                                    ### Function for calculating p_i 
                                    calculate_p <- function(my_data){
                                      my_data %>% 
                                        dplyr::mutate(p = ifelse(Y == 1,(Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data))/ Pr_Y1_theta(my_data), 
                                                                 ((1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data))/Pr_Y0_theta(my_data))) -> my_p_table
                                      return(my_p_table)
                                    }
                                    
                                    ### Function for calculating the expected complete log-likelihood  
                                    e_step_LL <- function(my_data){
                                      Y1 <- my_data %>% dplyr::filter(Y == 1)
                                      Y0 <- my_data %>% dplyr::filter(Y == 0)
                                      my_data %>% 
                                        dplyr::select(beta0,beta1) %>% 
                                        dplyr::distinct() %>% 
                                        as.numeric() -> my_beta
                                      Y1_X1 <- cbind(1, as.matrix(Y1['X1']))
                                      Y0_X1 <- cbind(1, as.matrix(Y0['X1']))
                                      
                                      LL_Y1 <- sum((Y1$p*(Y1$Y*Y1$fMi - log(1 + exp(Y1$fMi)))) + 
                                                     ((1-Y1$p)*(Y1$Y*Y1$epsilon - log(1+exp(Y1$epsilon)))) + 
                                                     (Y1$p*(Y1_X1 %*% my_beta) - log(1+exp(Y1_X1 %*% my_beta)))) 
                                      
                                      LL_Y0 <- sum((Y0$p)*(Y0$Y*Y0$fMi - log(1 + exp(Y0$fMi))) + 
                                                     ((1-Y0$p)*(Y0$Y*Y0$epsilon - log(1+exp(Y0$epsilon)))) + 
                                                     (Y0$p*Y0_X1 %*% my_beta - log(1+exp(Y0_X1 %*% my_beta)))) 
                                      
                                      Expected_LL = LL_Y1 + LL_Y0   
                                      return(Expected_LL)
                                    }
                                    
                                    ### Function for calculating the incomplete log-likelihood 
                                    incomplete_LL <- function(par, my_data){
                                      e = 0 # e is parameter argument 1 
                                      vBeta = c(par[1],par[2]) # my beta vector is parameter arguments 2 and 3 
                                      Y0 <- my_data %>% dplyr::filter(Y == 0)
                                      Y1 <- my_data %>% dplyr::filter(Y == 1)
                                      
                                      Mi0 = as.matrix(Y0['fMi']) # f_tilde for Y=0
                                      # Mi0 = log(Mi0)
                                      # add an intercept to the predictor variables
                                      #Mi0 = cbind(1, Mi0)
                                      mX0 = as.matrix(Y0['X1'])
                                      mX0 = cbind(1, mX0)
                                      
                                      Mi1 = as.matrix(Y1['fMi']) # f_tilde for Y = 1
                                      # Mi1 = log(Mi1) # modeling the log coverage
                                      #Mi1 = cbind(1, Mi1)
                                      mX1 = as.matrix(Y1['X1'])
                                      mX1 = cbind(1, mX1)
                                      
                                      LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                      (1-(exp(Mi0)/(1+exp(Mi0))))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                        sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                  (1/(1+exp(-(mX1 %*% vBeta))))*(exp(Mi1)/(1+exp(Mi1)))))
                                      return(LL)
                                    }
                                    
                                    ### Function for Maximization of vBeta
                                    vBeta <- function(my_data){
                                      ### logistic regression:
                                      # my_updated_beta <- glm(p~X1, data = my_data, family = binomial)
                                      ### Firth-penalized logistic regression:
                                      my_updated_beta <- logistf::logistf(p~X1, data = my_data)
                                      my_vBeta <-c(my_updated_beta$coefficients[[1]],my_updated_beta$coefficients[[2]])
                                      return(my_vBeta)
                                    }
                                    
                                    ### Function for maximization of vBeta Null model 
                                    vBeta_null <- function(my_data){
                                      ### logistic regression:
                                      # my_updated_beta <- glm(p~X1, data = my_data, family = binomial)
                                      ### Firth-penalized logistic regression:
                                      my_updated_beta <- logistf::logistf(p~1, data = my_data)
                                      my_vBeta <- my_updated_beta$coefficients[[1]]
                                      return(my_vBeta)
                                    }
                                    
                                    ### Function for Maximization of fMi 
                                    vfMi <- function(my_data){
                                      my_data %>% dplyr::arrange(coverage) -> my_data
                                      my_data$Y -> y 
                                      my_data$p -> w1
                                      n <- nrow(my_data)
                                      
                                      Atot <- cbind(1:(n-1),2:n) # define monotonicity 
                                      
                                      fit <- isotone::activeSet(isomat = Atot,
                                                                mySolver = fSolver,
                                                                fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/50)^2)),
                                                                gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (x/25)*sinh((x/50)^2),
                                                                y = y,
                                                                weights = w1)
                                      # a =  25
                                      # fit <- isotone::activeSet(isomat = Atot,
                                      #                           mySolver = fSolver,
                                      #                           fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/a)^4)),
                                      #                           gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (4*(x^3)/a)*sinh((x/50)^4),
                                      #                           y = y,
                                      #                           weights = w1)
                                      fit2 <- cbind(fit$y,fit$x,my_data$coverage) %>% as.data.frame()
                                      colnames(fit2) <- c("Y","fMi","coverage")
                                      return(fit2)
                                    }
                                    
                                    # if the Y vector has all 1's or 0's then pvalue = 1
                                    if(length(draw$Y %>% unique()) == 1){
                                      my_results <- c(0,0,0,0,0)
                                      my_pvalue <- 1
                                      #my_LRT <- 0
                                      
                                    } else if(length(draw$X1 %>% unique()) == 1) { 
                                      # if my covariate is all 1's or 0's then pvalue = 1 (consider something else to flag this data structure?)
                                      my_results <- c(0,0,0,0,0)
                                      my_pvalue <- 1
                                      #my_LRT <- 0
                                      
                                    } else {
                                      draw$epsilon <- 0
                                      num <- nrow(draw)
                                      
                                      ###### Need to estimate the betas for the full dataset first 
                                      
                                      init_list <- list()
                                      init_list[[1]] <- rep(1,num)
                                      #init_list[[2]] <- seq(-1,1,length = num)
                                      
                                      final_LL_list <- list() 
                                      
                                      for (k in 1:length(init_list)){  
                                        
                                        draw %>% dplyr::arrange(coverage) -> draw
                                        draw$fMi <- init_list[[k]] # Start with a uniform fMi's
                                        # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                        
                                        #t <- 100 # how many iterations for algorithm to converge
                                        my_data1 <- calculate_p(draw) # initial dataset calculation of p's 
                                        my_data1 <- my_data1 %>% 
                                          dplyr::mutate(iteration = 1)
                                        
                                        LL_list <- list() # List to store LL's 
                                        
                                        mydata_list <- list() # List to store datasets corresponding to LL's 
                                        
                                        ftilde_list <- list() 
                                        
                                        i <- 1 
                                        beta1_change <- 1
                                        beta0_change <- 1
                                        fmi_change <- 1
                                        rel_tol <- 1
                                        #abs_tol <- 1
                                        eval(as.symbol(paste0("my_data",i))) %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric() -> my_init_betas
                                        
                                        LL_list[i] <- incomplete_LL(my_init_betas,eval(as.symbol(paste0("my_data",i))))
                                        mydata_list[i] <- list(eval(as.symbol(paste0("my_data",i))))
                                        ftilde_list[i] <- -1 
                                        while (beta1_change > 0.0001 | beta0_change > 0.0001 | fmi_change > 0.01) {
                                          j <- i+1 # increment by i + 1 
                                          # I need to calculate the p's which I did outside this loop for the initial dataset
                                          # Then calculate my LL for that iteration 
                                          
                                          # Using p's I can update my parameters for the next dataset j to use 
                                          assign(paste0("nextbetas_",j), vBeta(eval(as.symbol(paste0("my_data",i)))))
                                          assign(paste0("nextfMis_",j), vfMi(eval(as.symbol(paste0("my_data",i)))))
                                          
                                          # Then I need to Drop the previous variables that will be replaced with new ones
                                          drop.cols <- c('fMi','beta0','beta1','p')
                                          assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",i))) %>% dplyr::select(-drop.cols))
                                          # Then Add the new betas and fMis to the dataframe 
                                          assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",j))) %>% 
                                                   dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetas_",j)))[1], 
                                                                 beta1 = eval(as.symbol(paste0("nextbetas_",j)))[2], 
                                                                 iteration = j))
                                          mergeCols <- c("coverage","Y")
                                          assign(paste0("my_data",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMis_",j))), eval(as.symbol(paste0("my_data",j))), by = mergeCols))
                                          #Now I have an updated dataframe of fMis and betas 
                                          
                                          #Next I need to Calculate the new P's for this 
                                          
                                          assign(paste0("my_data",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_data",j)))))
                                          
                                          mydata_list[j] <- list(eval(as.symbol(paste0("my_data",j))))
                                          # Then calculate my LL for that iteration 
                                          eval(as.symbol(paste0("my_data",j))) %>% 
                                            dplyr::select(beta0,beta1) %>% 
                                            dplyr::distinct() %>% 
                                            as.numeric() -> my_j_betas
                                          
                                          LL_list[j] <- incomplete_LL(my_j_betas,eval(as.symbol(paste0("my_data",j))))
                                          
                                          #  LL_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                          
                                          #abs_tol <- abs(LL_list[[j]] - LL_list[[i]])
                                          #rel_tol <- (abs(LL_list[[j]] - LL_list[[i]]))/abs(LL_list[[i]])
                                          
                                          betas_next <- mydata_list[[j]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          betas_previous <- mydata_list[[i]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          
                                          beta0_change <- abs(betas_next[[1]]-betas_previous[[1]])
                                          beta1_change <- abs(betas_next[[2]]-betas_previous[[2]])
                                          
                                          fmis_next <- mydata_list[[j]] %>% dplyr::select(fMi) 
                                          fmis_previous <- mydata_list[[i]] %>% dplyr::select(fMi)
                                          
                                          fmis_intermediate <- abs(fmis_next-fmis_previous)/abs(fmis_previous)
                                          max(fmis_intermediate)
                                          fmi_change <- max(fmis_intermediate)
                                          
                                          ftilde_list[j] <- fmi_change
                                          
                                          i = i + 1
                                          
                                          print(j)
                                          if (j == 500){
                                            rel_tol <- -1 
                                            beta1_change <- -1
                                            beta0_change <- -1
                                            fmi_change <- -1 
                                            #LL_list[[j]] <- "max iterations reached"
                                            # do.call("rbind", ftilde_list) -> ftilde_table
                                            # ftilde_table %>% 
                                            #   dplyr::as_data_frame() %>% 
                                            #   dplyr::mutate(iteration = 1:dplyr::n()) %>% 
                                            #   dplyr::filter(iteration > 1) -> ftilde_change_plot
                                            # ftilde_change_plot[which.min(ftilde_change_plot$V1),][[1]] -> ftilde_list[[j]]
                                          }
                                        }
                                        length(mydata_list) -> em_iteration
                                        
                                        mydata_list[[em_iteration]] %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric() -> my_em_betas
                                        
                                        cbind(my_em_betas[[1]], my_em_betas[[2]],LL_list[[em_iteration]], 
                                              em_iteration[[1]],ftilde_list[[em_iteration]]) -> final_LL_list[[k]]
                                      }
                                      
                                      do.call(rbind, final_LL_list) %>% 
                                        dplyr::as_tibble() %>% 
                                        dplyr::filter(!is.na(V1)) %>% 
                                        dplyr::distinct() -> results_tibble
                                      results_tibble[which.min(results_tibble$V5),] -> my_results 
                                      if (dim(my_results)[[1]] == 0){
                                        print("unrestricted model didn't converge")
                                        my_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                      } else {
                                        print("unrestricted model completed!")
                                      }
                                      
                                      
                                      
                                      
                                      #### Here starts the Wasserman related code ######
                                      
                                      picked <- sample(seq_len(nrow(draw)), size = num/2, replace = FALSE)
                                      
                                      D_1 <- draw[picked,]
                                      D_0 <- draw[-picked,]
                                      
                                      init_list <- list()
                                      init_list[[1]] <- rep(1,nrow(D_1))
                                      #init_list[[2]] <- seq(-1,1,length = nrow(D_1))
                                      
                                      final_LL_list_D1 <- list() 
                                      
                                      for (k in 1:length(init_list)){  
                                        
                                        D_1 %>% dplyr::arrange(coverage) -> D_1
                                        D_1$fMi <- init_list[[k]] # Start with a uniform fMi's
                                        # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                        
                                        #t <- 100 # how many iterations for algorithm to converge
                                        my_data1 <- calculate_p(D_1) # initial dataset calculation of p's 
                                        my_data1 <- my_data1 %>% 
                                          dplyr::mutate(iteration = 1)
                                        
                                        LL_list_theta_D1_fullmodel <- list() # List to store LL's 
                                        
                                        mydata_list_D1 <- list() # List to store datasets corresponding to LL's 
                                        
                                        ftilde_list_D1 <- list() 
                                        
                                        i <- 1 
                                        beta1_change <- 1
                                        beta0_change <- 1
                                        fmi_change <- 1
                                        rel_tol <- 1
                                        #abs_tol <- 1
                                        eval(as.symbol(paste0("my_data",i))) %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric() -> my_init_betas
                                        
                                        LL_list_theta_D1_fullmodel[i] <- incomplete_LL(my_init_betas,D_0)
                                        mydata_list_D1[i] <- list(eval(as.symbol(paste0("my_data",i))))
                                        ftilde_list_D1[i] <- -1 
                                        
                                        while (beta1_change > 0.0001 | beta0_change > 0.0001 | fmi_change > 0.01) {
                                          j <- i+1 # increment by i + 1 
                                          # I need to calculate the p's which I did outside this loop for the initial dataset
                                          # Then calculate my LL for that iteration 
                                          
                                          # Using p's I can update my parameters for the next dataset j to use 
                                          assign(paste0("nextbetas_",j), vBeta(eval(as.symbol(paste0("my_data",i)))))
                                          assign(paste0("nextfMis_",j), vfMi(eval(as.symbol(paste0("my_data",i)))))
                                          
                                          # Then I need to Drop the previous variables that will be replaced with new ones
                                          drop.cols <- c('fMi','beta0','beta1','p')
                                          assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",i))) %>% dplyr::select(-drop.cols))
                                          # Then Add the new betas and fMis to the dataframe 
                                          assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",j))) %>% 
                                                   dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetas_",j)))[1], 
                                                                 beta1 = eval(as.symbol(paste0("nextbetas_",j)))[2], 
                                                                 iteration = j))
                                          mergeCols <- c("coverage","Y")
                                          assign(paste0("my_data",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMis_",j))), eval(as.symbol(paste0("my_data",j))), by = mergeCols))
                                          #Now I have an updated dataframe of fMis and betas 
                                          
                                          #Next I need to Calculate the new P's for this 
                                          
                                          assign(paste0("my_data",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_data",j)))))
                                          
                                          mydata_list_D1[j] <- list(eval(as.symbol(paste0("my_data",j))))
                                          # Then calculate my LL for that iteration 
                                          eval(as.symbol(paste0("my_data",j))) %>% 
                                            dplyr::select(beta0,beta1) %>% 
                                            dplyr::distinct() %>% 
                                            as.numeric() -> my_j_betas
                                          
                                          LL_list_theta_D1_fullmodel[j] <- incomplete_LL(my_j_betas,D_0)
                                          
                                          #  LL_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                          
                                          #abs_tol <- abs(LL_list[[j]] - LL_list[[i]])
                                          #rel_tol <- (abs(LL_list[[j]] - LL_list[[i]]))/abs(LL_list[[i]])
                                          
                                          betas_next <- mydata_list_D1[[j]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          betas_previous <- mydata_list_D1[[i]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          
                                          beta0_change <- abs(betas_next[[1]]-betas_previous[[1]])
                                          beta1_change <- abs(betas_next[[2]]-betas_previous[[2]])
                                          
                                          fmis_next <- mydata_list_D1[[j]] %>% dplyr::select(fMi) 
                                          fmis_previous <- mydata_list_D1[[i]] %>% dplyr::select(fMi)
                                          
                                          fmis_intermediate <- abs(fmis_next-fmis_previous)/abs(fmis_previous)
                                          max(fmis_intermediate)
                                          fmi_change <- max(fmis_intermediate)
                                          
                                          ftilde_list_D1[j] <- fmi_change
                                          
                                          i = i + 1
                                          
                                          print(j)
                                          if (j == 500){
                                            rel_tol <- -1 
                                            beta1_change <- -1
                                            beta0_change <- -1
                                            fmi_change <- -1 
                                            #LL_list[[j]] <- "max iterations reached"
                                            # do.call("rbind", ftilde_list) -> ftilde_table
                                            # ftilde_table %>% 
                                            #   dplyr::as_data_frame() %>% 
                                            #   dplyr::mutate(iteration = 1:dplyr::n()) %>% 
                                            #   dplyr::filter(iteration > 1) -> ftilde_change_plot
                                            # ftilde_change_plot[which.min(ftilde_change_plot$V1),][[1]] -> ftilde_list[[j]]
                                          }
                                        }
                                        length(mydata_list_D1) -> em_iteration
                                        
                                        mydata_list_D1[[em_iteration]] %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric() -> my_em_betas
                                        
                                        cbind(my_em_betas[[1]], my_em_betas[[2]],LL_list_theta_D1_fullmodel[[em_iteration]], 
                                              em_iteration[[1]],ftilde_list_D1[[em_iteration]]) -> final_LL_list_D1[[k]]
                                      }
                                      
                                      do.call(rbind, final_LL_list_D1) %>% 
                                        dplyr::as_tibble() %>% 
                                        dplyr::filter(!is.na(V1)) %>% 
                                        dplyr::distinct() -> results_tibble_D1
                                      results_tibble_D1[which.min(results_tibble_D1$V5),] -> my_results_D1 
                                      if (dim(my_results_D1)[[1]] == 0){
                                        print("unrestricted model didn't converge")
                                        my_results_D1 <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                      } else {
                                        print("unrestricted model completed!")
                                      }
                                      
                                      
                                      ############# for the restricted/null model ###############
                                      
                                      init_list <- list()
                                      init_list[[1]] <- rep(1,nrow(D_1))
                                      #init_list[[2]] <- seq(-1,1,length = nrow(D_1))
                                      
                                      
                                      
                                      LL_list_theta_D1_restricted <- list() 
                                      
                                      for (k in 1:length(init_list)){  
                                        
                                        D_1 %>% dplyr::arrange(coverage) -> D_1
                                        D_1$fMi <- init_list[[k]] # Start with a uniform fMi's
                                        # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                        D_1 %>% 
                                          dplyr::select(-beta1) %>% 
                                          dplyr::mutate(beta1 = 0) -> null_draw
                                        #t <- 100 # how many iterations for algorithm to converge
                                        my_nulldata1 <- calculate_p(null_draw) # initial dataset calculation of p's 
                                        my_nulldata1 <- my_nulldata1 %>% 
                                          dplyr::mutate(iteration = 1)
                                        
                                        LL_restricted_list <- list() # List to store LL's 
                                        
                                        mydata_restricted_list <- list() # List to store datasets corresponding to LL's 
                                        
                                        ftilde_restricted_list <- list() 
                                        
                                        i <- 1 
                                        beta0_null_change <- 1
                                        fmi_null_change <- 1
                                        #rel_null_tol <- 1
                                        #abs_tol <- 1
                                        eval(as.symbol(paste0("my_nulldata",i))) %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric()  -> my_init_betas
                                        
                                        LL_list_theta_D1_restricted[i] <- incomplete_LL(my_init_betas,D_0)
                                        mydata_restricted_list[i] <- list(eval(as.symbol(paste0("my_nulldata",i))))
                                        ftilde_restricted_list[i] <- -1 
                                        
                                        while (beta0_null_change > 0.0001 | fmi_null_change > 0.01) {
                                          j <- i+1 # increment by i + 1 
                                          # I need to calculate the p's which I did outside this loop for the initial dataset
                                          # Then calculate my LL for that iteration 
                                          
                                          # Using p's I can update my parameters for the next dataset j to use 
                                          assign(paste0("nextbetasnull_",j), vBeta_null(eval(as.symbol(paste0("my_nulldata",i)))))
                                          assign(paste0("nextfMisnull_",j), vfMi(eval(as.symbol(paste0("my_nulldata",i)))))
                                          
                                          # Then I need to Drop the previous variables that will be replaced with new ones
                                          drop.cols <- c('fMi','beta0','beta1','p')
                                          assign(paste0("my_nulldata",j), eval(as.symbol(paste0("my_nulldata",i))) %>% dplyr::select(-drop.cols))
                                          # Then Add the new betas and fMis to the dataframe 
                                          assign(paste0("my_nulldata",j), eval(as.symbol(paste0("my_nulldata",j))) %>% 
                                                   dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetasnull_",j)))[1],
                                                                 beta1 = 0 , 
                                                                 iteration = j))
                                          mergeCols <- c("coverage","Y")
                                          assign(paste0("my_nulldata",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMisnull_",j))), eval(as.symbol(paste0("my_nulldata",j))), by = mergeCols))
                                          #Now I have an updated dataframe of fMis and betas 
                                          
                                          #Next I need to Calculate the new P's for this 
                                          
                                          assign(paste0("my_nulldata",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_nulldata",j)))))
                                          
                                          mydata_restricted_list[j] <- list(eval(as.symbol(paste0("my_nulldata",j))))
                                          # Then calculate my LL for that iteration 
                                          eval(as.symbol(paste0("my_nulldata",j))) %>% 
                                            dplyr::select(beta0,beta1) %>% 
                                            dplyr::distinct() %>% 
                                            as.numeric() -> my_j_betas 
                                          
                                          LL_list_theta_D1_restricted[j] <- incomplete_LL(my_j_betas,D_0)
                                          
                                          #LL_restricted_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                          
                                          #abs_tol <- abs(LL_restricted_list[[j]] - LL_restricted_list[[i]])
                                          #rel_tol <- (abs(LL_restricted_list[[j]] - LL_restricted_list[[i]]))/abs(LL_restricted_list[[i]])
                                          
                                          betas_null_next <- mydata_restricted_list[[j]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          betas_null_previous <- mydata_restricted_list[[i]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          
                                          beta0_null_change <- abs(betas_null_next[[1]]-betas_null_previous[[1]])
                                          beta1_null_change <- abs(betas_null_next[[2]]-betas_null_previous[[2]])
                                          
                                          fmis_null_next <- mydata_restricted_list[[j]] %>% dplyr::select(fMi) 
                                          fmis_null_previous <- mydata_restricted_list[[i]] %>% dplyr::select(fMi)
                                          
                                          fmis_null_intermediate <- abs(fmis_null_next-fmis_null_previous)/abs(fmis_null_previous)
                                          max(fmis_null_intermediate)
                                          fmi_null_change <- max(fmis_null_intermediate)
                                          
                                          ftilde_restricted_list[j] <- fmi_null_change
                                          
                                          i = i + 1
                                          
                                          print(j)
                                          if (j == 500){
                                            rel_tol <- -1 
                                            beta1_null_change <- -1
                                            beta0_null_change <- -1
                                            fmi_null_change <- -1 
                                            #LL_list[[j]] <- "max iterations reached"
                                            # do.call("rbind", ftilde_restricted_list) -> ftilde_table
                                            # ftilde_table %>% 
                                            #   dplyr::as_data_frame() %>% 
                                            #   dplyr::mutate(iteration = 1:dplyr::n()) %>% 
                                            #   dplyr::filter(iteration > 1) -> ftilde_change_plot
                                            # ftilde_change_plot[which.min(ftilde_change_plot$V1),][[1]] -> ftilde_list[[j]]
                                          }
                                        }
                                        
                                        length(mydata_restricted_list) -> em_restricted_iteration
                                        
                                        mydata_restricted_list[[em_restricted_iteration]] %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric() -> my_em_restricted_betas
                                        
                                        cbind(my_em_restricted_betas[[1]], my_em_restricted_betas[[2]],LL_list_theta_D1_restricted[[em_restricted_iteration]], 
                                              em_restricted_iteration[[1]],ftilde_restricted_list[[em_restricted_iteration]]) -> LL_restricted_list[[k]]
                                      }
                                      
                                      do.call(rbind, LL_restricted_list) %>% 
                                        dplyr::as_tibble() %>% 
                                        dplyr::filter(!is.na(V1)) %>% 
                                        dplyr::distinct() -> restricted_results_tibble
                                      restricted_results_tibble[which.min(results_tibble$V5),] -> my_restricted_results 
                                      if (dim(my_restricted_results)[[1]] == 0){
                                        print("restricted model didn't converge")
                                        my_restricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                      } else {
                                        print("restricted model completed!")
                                      }
                                      
                                      U_n = exp(my_results_D1$V3 - my_restricted_results$V3)
                                      
                                      
                                      #########################################################
                                      ##### Calculating U_n_swap ##############################
                                      #########################################################
                                      # D_1 
                                      # D_0 
                                      
                                      init_list_swap <- list()
                                      init_list_swap[[1]] <- rep(1,nrow(D_0))
                                      #init_list[[2]] <- seq(-1,1,length = nrow(D_1))
                                      
                                      final_LL_list_swap <- list() 
                                      
                                      for (k in 1:length(init_list_swap)){  
                                        
                                        D_0 %>% dplyr::arrange(coverage) -> D_0
                                        D_0$fMi <- init_list_swap[[k]] # Start with a uniform fMi's
                                        # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                        
                                        #t <- 100 # how many iterations for algorithm to converge
                                        my_data1 <- calculate_p(D_0) # initial dataset calculation of p's 
                                        my_data1 <- my_data1 %>% 
                                          dplyr::mutate(iteration = 1)
                                        
                                        LL_list_theta_D0_fullmodel <- list() # List to store LL's 
                                        
                                        mydata_list_swap <- list() # List to store datasets corresponding to LL's 
                                        
                                        ftilde_list_swap <- list() 
                                        
                                        i <- 1 
                                        beta1_change <- 1
                                        beta0_change <- 1
                                        fmi_change <- 1
                                        rel_tol <- 1
                                        #abs_tol <- 1
                                        eval(as.symbol(paste0("my_data",i))) %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric() -> my_init_betas
                                        
                                        LL_list_theta_D0_fullmodel[i] <- incomplete_LL(my_init_betas,D_1)
                                        mydata_list_swap[i] <- list(eval(as.symbol(paste0("my_data",i))))
                                        ftilde_list_swap[i] <- -1 
                                        
                                        while (beta1_change > 0.0001 | beta0_change > 0.0001 | fmi_change > 0.01) {
                                          j <- i+1 # increment by i + 1 
                                          # I need to calculate the p's which I did outside this loop for the initial dataset
                                          # Then calculate my LL for that iteration 
                                          
                                          # Using p's I can update my parameters for the next dataset j to use 
                                          assign(paste0("nextbetas_",j), vBeta(eval(as.symbol(paste0("my_data",i)))))
                                          assign(paste0("nextfMis_",j), vfMi(eval(as.symbol(paste0("my_data",i)))))
                                          
                                          # Then I need to Drop the previous variables that will be replaced with new ones
                                          drop.cols <- c('fMi','beta0','beta1','p')
                                          assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",i))) %>% dplyr::select(-drop.cols))
                                          # Then Add the new betas and fMis to the dataframe 
                                          assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",j))) %>% 
                                                   dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetas_",j)))[1], 
                                                                 beta1 = eval(as.symbol(paste0("nextbetas_",j)))[2], 
                                                                 iteration = j))
                                          mergeCols <- c("coverage","Y")
                                          assign(paste0("my_data",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMis_",j))), eval(as.symbol(paste0("my_data",j))), by = mergeCols))
                                          #Now I have an updated dataframe of fMis and betas 
                                          
                                          #Next I need to Calculate the new P's for this 
                                          
                                          assign(paste0("my_data",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_data",j)))))
                                          
                                          mydata_list_swap[j] <- list(eval(as.symbol(paste0("my_data",j))))
                                          # Then calculate my LL for that iteration 
                                          eval(as.symbol(paste0("my_data",j))) %>% 
                                            dplyr::select(beta0,beta1) %>% 
                                            dplyr::distinct() %>% 
                                            as.numeric() -> my_j_betas
                                          
                                          LL_list_theta_D0_fullmodel[j] <- incomplete_LL(my_j_betas,D_1)
                                          
                                          #  LL_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                          
                                          #abs_tol <- abs(LL_list[[j]] - LL_list[[i]])
                                          #rel_tol <- (abs(LL_list[[j]] - LL_list[[i]]))/abs(LL_list[[i]])
                                          
                                          betas_next <- mydata_list_swap[[j]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          betas_previous <- mydata_list_swap[[i]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          
                                          beta0_change <- abs(betas_next[[1]]-betas_previous[[1]])
                                          beta1_change <- abs(betas_next[[2]]-betas_previous[[2]])
                                          
                                          fmis_next <- mydata_list_swap[[j]] %>% dplyr::select(fMi) 
                                          fmis_previous <- mydata_list_swap[[i]] %>% dplyr::select(fMi)
                                          
                                          fmis_intermediate <- abs(fmis_next-fmis_previous)/abs(fmis_previous)
                                          max(fmis_intermediate)
                                          fmi_change <- max(fmis_intermediate)
                                          
                                          ftilde_list_swap[j] <- fmi_change
                                          
                                          i = i + 1
                                          
                                          print(j)
                                          if (j == 500){
                                            rel_tol <- -1 
                                            beta1_change <- -1
                                            beta0_change <- -1
                                            fmi_change <- -1 
                                            #LL_list[[j]] <- "max iterations reached"
                                            # do.call("rbind", ftilde_list) -> ftilde_table
                                            # ftilde_table %>% 
                                            #   dplyr::as_data_frame() %>% 
                                            #   dplyr::mutate(iteration = 1:dplyr::n()) %>% 
                                            #   dplyr::filter(iteration > 1) -> ftilde_change_plot
                                            # ftilde_change_plot[which.min(ftilde_change_plot$V1),][[1]] -> ftilde_list[[j]]
                                          }
                                        }
                                        length(mydata_list_swap) -> em_iteration
                                        
                                        mydata_list_swap[[em_iteration]] %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric() -> my_em_betas
                                        
                                        cbind(my_em_betas[[1]], my_em_betas[[2]],LL_list_theta_D0_fullmodel[[em_iteration]], 
                                              em_iteration[[1]],ftilde_list_swap[[em_iteration]]) -> final_LL_list_swap[[k]]
                                      }
                                      
                                      do.call(rbind, final_LL_list_swap) %>% 
                                        dplyr::as_tibble() %>% 
                                        dplyr::filter(!is.na(V1)) %>% 
                                        dplyr::distinct() -> results_tibble_swap
                                      results_tibble_swap[which.min(results_tibble_swap$V5),] -> my_results_swap 
                                      if (dim(my_results_swap)[[1]] == 0){
                                        print("unrestricted model swap didn't converge")
                                        my_results_swap <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                      } else {
                                        print("unrestricted model swap completed!")
                                      }
                                      
                                      
                                      ############# for the restricted/null model ###############
                                      
                                      init_list_swap <- list()
                                      init_list_swap[[1]] <- rep(1,nrow(D_0))
                                      #init_list[[2]] <- seq(-1,1,length = nrow(D_1))
                                      
                                      
                                      LL_restricted_list_swap <- list()
                                      
                                      
                                      for (k in 1:length(init_list_swap)){  
                                        
                                        D_0 %>% dplyr::arrange(coverage) -> D_0
                                        D_0$fMi <- init_list_swap[[k]] # Start with a uniform fMi's
                                        # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                        D_0 %>% 
                                          dplyr::select(-beta1) %>% 
                                          dplyr::mutate(beta1 = 0) -> null_draw
                                        #t <- 100 # how many iterations for algorithm to converge
                                        my_nulldata1 <- calculate_p(null_draw) # initial dataset calculation of p's 
                                        my_nulldata1 <- my_nulldata1 %>% 
                                          dplyr::mutate(iteration = 1)
                                        
                                        LL_list_theta_D0_restricted <- list() # List to store all LL's 
                                        
                                        mydata_restricted_list_swap <- list() # List to store datasets corresponding to LL's 
                                        
                                        ftilde_restricted_list_swap <- list() 
                                        
                                        i <- 1 
                                        beta0_null_change <- 1
                                        fmi_null_change <- 1
                                        #rel_null_tol <- 1
                                        #abs_tol <- 1
                                        eval(as.symbol(paste0("my_nulldata",i))) %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric()  -> my_init_betas
                                        
                                        LL_list_theta_D0_restricted[i] <- incomplete_LL(my_init_betas,D_1)
                                        mydata_restricted_list_swap[i] <- list(eval(as.symbol(paste0("my_nulldata",i))))
                                        ftilde_restricted_list_swap[i] <- -1 
                                        
                                        while (beta0_null_change > 0.0001 | fmi_null_change > 0.01) {
                                          j <- i+1 # increment by i + 1 
                                          # I need to calculate the p's which I did outside this loop for the initial dataset
                                          # Then calculate my LL for that iteration 
                                          
                                          # Using p's I can update my parameters for the next dataset j to use 
                                          assign(paste0("nextbetasnull_",j), vBeta_null(eval(as.symbol(paste0("my_nulldata",i)))))
                                          assign(paste0("nextfMisnull_",j), vfMi(eval(as.symbol(paste0("my_nulldata",i)))))
                                          
                                          # Then I need to Drop the previous variables that will be replaced with new ones
                                          drop.cols <- c('fMi','beta0','beta1','p')
                                          assign(paste0("my_nulldata",j), eval(as.symbol(paste0("my_nulldata",i))) %>% dplyr::select(-drop.cols))
                                          # Then Add the new betas and fMis to the dataframe 
                                          assign(paste0("my_nulldata",j), eval(as.symbol(paste0("my_nulldata",j))) %>% 
                                                   dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetasnull_",j)))[1],
                                                                 beta1 = 0 , 
                                                                 iteration = j))
                                          mergeCols <- c("coverage","Y")
                                          assign(paste0("my_nulldata",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMisnull_",j))), eval(as.symbol(paste0("my_nulldata",j))), by = mergeCols))
                                          #Now I have an updated dataframe of fMis and betas 
                                          
                                          #Next I need to Calculate the new P's for this 
                                          
                                          assign(paste0("my_nulldata",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_nulldata",j)))))
                                          
                                          mydata_restricted_list_swap[j] <- list(eval(as.symbol(paste0("my_nulldata",j))))
                                          # Then calculate my LL for that iteration 
                                          eval(as.symbol(paste0("my_nulldata",j))) %>% 
                                            dplyr::select(beta0,beta1) %>% 
                                            dplyr::distinct() %>% 
                                            as.numeric() -> my_j_betas 
                                          
                                          LL_list_theta_D0_restricted[j] <- incomplete_LL(my_j_betas,D_1)
                                          
                                          #LL_restricted_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                          
                                          #abs_tol <- abs(LL_restricted_list[[j]] - LL_restricted_list[[i]])
                                          #rel_tol <- (abs(LL_restricted_list[[j]] - LL_restricted_list[[i]]))/abs(LL_restricted_list[[i]])
                                          
                                          betas_null_next <- mydata_restricted_list_swap[[j]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          betas_null_previous <- mydata_restricted_list_swap[[i]] %>% dplyr::select(beta0,beta1) %>% dplyr::distinct() %>% as.numeric()
                                          
                                          beta0_null_change <- abs(betas_null_next[[1]]-betas_null_previous[[1]])
                                          beta1_null_change <- abs(betas_null_next[[2]]-betas_null_previous[[2]])
                                          
                                          fmis_null_next <- mydata_restricted_list_swap[[j]] %>% dplyr::select(fMi) 
                                          fmis_null_previous <- mydata_restricted_list_swap[[i]] %>% dplyr::select(fMi)
                                          
                                          fmis_null_intermediate <- abs(fmis_null_next-fmis_null_previous)/abs(fmis_null_previous)
                                          max(fmis_null_intermediate)
                                          fmi_null_change <- max(fmis_null_intermediate)
                                          
                                          ftilde_restricted_list_swap[j] <- fmi_null_change
                                          
                                          i = i + 1
                                          
                                          print(j)
                                          if (j == 500){
                                            rel_tol <- -1 
                                            beta1_null_change <- -1
                                            beta0_null_change <- -1
                                            fmi_null_change <- -1 
                                            #LL_list[[j]] <- "max iterations reached"
                                            # do.call("rbind", ftilde_restricted_list) -> ftilde_table
                                            # ftilde_table %>% 
                                            #   dplyr::as_data_frame() %>% 
                                            #   dplyr::mutate(iteration = 1:dplyr::n()) %>% 
                                            #   dplyr::filter(iteration > 1) -> ftilde_change_plot
                                            # ftilde_change_plot[which.min(ftilde_change_plot$V1),][[1]] -> ftilde_list[[j]]
                                          }
                                        }
                                        
                                        length(mydata_restricted_list_swap) -> em_restricted_iteration
                                        
                                        mydata_restricted_list_swap[[em_restricted_iteration]] %>% 
                                          dplyr::select(beta0,beta1) %>% 
                                          dplyr::distinct() %>% 
                                          as.numeric() -> my_em_restricted_betas
                                        
                                        cbind(my_em_restricted_betas[[1]], my_em_restricted_betas[[2]],LL_list_theta_D0_restricted[[em_restricted_iteration]], 
                                              em_restricted_iteration[[1]],ftilde_restricted_list_swap[[em_restricted_iteration]]) -> LL_restricted_list_swap[[k]]
                                      }
                                      
                                      do.call(rbind, LL_restricted_list_swap) %>% 
                                        dplyr::as_tibble() %>% 
                                        dplyr::filter(!is.na(V1)) %>% 
                                        dplyr::distinct() -> restricted_results_tibble_swap
                                      restricted_results_tibble_swap[which.min(restricted_results_tibble_swap$V5),] -> my_restricted_results_swap
                                      if (dim(my_restricted_results_swap)[[1]] == 0){
                                        print("restricted model swap didn't converge")
                                        my_restricted_results_swap <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA)
                                      } else {
                                        print("restricted model swap completed!")
                                      }
                                      
                                      
                                      U_n_swap = exp(my_results_swap$V3 - my_restricted_results_swap$V3) 
                                      
                                      my_pvalue = 1/((U_n + U_n_swap)/2)
                                      
                                    }
                                    
                                    output <- list("beta0" = my_results[[1]],
                                                   "beta1" = my_results[[2]], 
                                                   "iterations" = my_results[[4]], 
                                                   "ftilde_change" = my_results[[5]], 
                                                   "pvalue"= my_pvalue)
                                  })


david_method_LLconvergence <- new_method("david_method_LLconvergence","david_method_LLconvergence",
                                         method = function(model,draw){
                                           ### Functions for the 3 hierarchical model components 
                                           # Probability of Y = 1 given Lambda = 1
                                           library(isotone)
                                           Pr_Y1_Lambda1 <- function(my_fMi) {
                                             my_Pr_Y1_Lambda1 <- exp(my_fMi)/(1+exp(my_fMi))
                                             return(my_Pr_Y1_Lambda1)
                                           }
                                           # Probability of Lambda = 1 
                                           Pr_Lambda1 <- function(my_data) {
                                             mX1 = as.matrix(my_data['X1'])
                                             mX1 = cbind(1, mX1)
                                             my_data %>% 
                                               dplyr::select(beta0,beta1) %>% 
                                               dplyr::distinct() %>% 
                                               as.numeric() -> my_beta
                                             my_Pr_Lambda1 <- exp(mX1 %*% my_beta)/(1+exp(mX1 %*% my_beta))
                                             return(my_Pr_Lambda1)
                                           }
                                           # Probability of Y = 1 given Lambda = 0 
                                           Pr_Y1_Lambda0 <- function(my_epsilon) {
                                             #my_Pr_Y1_Lambda0 <- exp(my_epsilon)/(1+exp(my_epsilon))
                                             # for now I'm going to set this parameterization to 0 
                                             my_Pr_Y1_Lambda0 <- 0 
                                             return(my_Pr_Y1_Lambda0)
                                           }
                                           
                                           ### Functions for calculating the denominators of p_i 
                                           # Prob of Y = 1 given parameters theta 
                                           Pr_Y1_theta <- function(my_data){
                                             my_Pr_Y1_theta <- Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data) + Pr_Y1_Lambda0(my_data$epsilon)*(1-Pr_Lambda1(my_data))
                                             return(my_Pr_Y1_theta)
                                           }
                                           #Prob of Y = 0 given parameters theta 
                                           Pr_Y0_theta <- function(my_data){
                                             my_Pr_Y0_theta <- (1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data) + (1-Pr_Y1_Lambda0(my_data$epsilon))*(1-Pr_Lambda1(my_data))
                                             return(my_Pr_Y0_theta)
                                           }
                                           
                                           ### Function for calculating p_i 
                                           calculate_p <- function(my_data){
                                             my_data %>% 
                                               dplyr::mutate(p = ifelse(Y == 1,(Pr_Y1_Lambda1(my_data$fMi)*Pr_Lambda1(my_data))/ Pr_Y1_theta(my_data), 
                                                                        ((1-Pr_Y1_Lambda1(my_data$fMi))*Pr_Lambda1(my_data))/Pr_Y0_theta(my_data))) -> my_p_table
                                             return(my_p_table)
                                           }
                                           
                                           ### Function for calculating the complete log-likelihood  
                                           e_step_LL <- function(my_data){
                                             Y1 <- my_data %>% dplyr::filter(Y == 1)
                                             Y0 <- my_data %>% dplyr::filter(Y == 0)
                                             my_data %>% 
                                               dplyr::select(beta0,beta1) %>% 
                                               dplyr::distinct() %>% 
                                               as.numeric() -> my_beta
                                             Y1_X1 <- cbind(1, as.matrix(Y1['X1']))
                                             Y0_X1 <- cbind(1, as.matrix(Y0['X1']))
                                             
                                             LL_Y1 <- sum((Y1$p*(Y1$Y*Y1$fMi - log(1 + exp(Y1$fMi)))) + 
                                                            ((1-Y1$p)*(Y1$Y*Y1$epsilon - log(1+exp(Y1$epsilon)))) + 
                                                            (Y1$p*(Y1_X1 %*% my_beta) - log(1+exp(Y1_X1 %*% my_beta)))) 
                                             
                                             LL_Y0 <- sum((Y0$p)*(Y0$Y*Y0$fMi - log(1 + exp(Y0$fMi))) + 
                                                            ((1-Y0$p)*(Y0$Y*Y0$epsilon - log(1+exp(Y0$epsilon)))) + 
                                                            (Y0$p*Y0_X1 %*% my_beta - log(1+exp(Y0_X1 %*% my_beta)))) 
                                             
                                             Expected_LL = LL_Y1 + LL_Y0   
                                             return(Expected_LL)
                                           }
                                           
                                           ### Function for Maximization of vBeta
                                           vBeta <- function(my_data){
                                             ### logistic regression:
                                             # my_updated_beta <- glm(p~X1, data = my_data, family = binomial)
                                             ### Firth-penalized logistic regression:
                                             my_updated_beta <- logistf::logistf(p~X1, data = my_data)
                                             my_vBeta <-c(my_updated_beta$coefficients[[1]],my_updated_beta$coefficients[[2]])
                                             return(my_vBeta)
                                           }
                                           
                                           ### Function for Maximization of fMi 
                                           vfMi <- function(my_data){
                                             my_data %>% dplyr::arrange(coverage) -> my_data
                                             my_data$Y -> y 
                                             my_data$p -> w1
                                             n <- nrow(my_data)
                                             
                                             Atot <- cbind(1:(n-1),2:n) # define monotonicity 
                                             
                                             fit <- isotone::activeSet(isomat = Atot,
                                                                       mySolver = fSolver,
                                                                       fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/50)^2)),
                                                                       gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (x/25)*sinh((x/50)^2),
                                                                       y = y,
                                                                       weights = w1)
                                             # a =  25
                                             # fit <- isotone::activeSet(isomat = Atot,
                                             #                           mySolver = fSolver,
                                             #                           fobj = function(x) -1*sum(w1*(y*x-log(1+exp(x)))) + sum(cosh((x/a)^4)),
                                             #                           gobj = function(x) -1*(w1*(y-(exp(x)/(1+exp(x))))) + (4*(x^3)/a)*sinh((x/50)^4),
                                             #                           y = y,
                                             #                           weights = w1)
                                             fit2 <- cbind(fit$y,fit$x,my_data$coverage) %>% as.data.frame()
                                             colnames(fit2) <- c("Y","fMi","coverage")
                                             return(fit2)
                                           }
                                           
                                           # Okay I need to create multiple starts that have LL's stored in a list and the set of results that has 
                                           # the highest LL will be the chosen one 
                                           num <- nrow(draw)
                                           draw$epsilon <- 0
                                           
                                           init_list <- list()
                                           init_list[[1]] <- rep(1,num)
                                           init_list[[2]] <- rep(0.1,num)
                                           init_list[[3]] <- rep(2,num)
                                           init_list[[4]] <- rep(-1,num)
                                           
                                           
                                           final_LL_list <- list() 
                                           
                                           for (k in 1:length(init_list)){  
                                             
                                             draw$fMi <- init_list[[k]] # Start with a uniform fMi's
                                             # my_data <- calculate_p(my_data) # initial dataset calculation of p's 
                                             
                                             #t <- 100 # how many iterations for algorithm to converge
                                             my_data1 <- calculate_p(draw) # initial dataset calculation of p's 
                                             my_data1 <- my_data1 %>% 
                                               dplyr::mutate(iteration = 1)
                                             LL_list <- list() # List to store LL's 
                                             
                                             mydata_list <- list() # List to store datasets corresponding to LL's 
                                             
                                             i <- 1 
                                             beta1_change <- 1
                                             beta0_change <- 1
                                             rel_tol <- 1
                                             #abs_tol <- 1
                                             LL_list[i] <- e_step_LL(eval(as.symbol(paste0("my_data",i))))
                                             mydata_list[i] <- list(eval(as.symbol(paste0("my_data",i))))
                                             
                                             while (rel_tol > 0.00001) {
                                               j <- i+1 # increment by i + 1 
                                               # I need to calculate the p's which I did outside this loop for the initial dataset
                                               # Then calculate my LL for that iteration 
                                               
                                               # Using p's I can update my parameters for the next dataset j to use 
                                               assign(paste0("nextbetas_",j), vBeta(eval(as.symbol(paste0("my_data",i)))))
                                               assign(paste0("nextfMis_",j), vfMi(eval(as.symbol(paste0("my_data",i)))))
                                               
                                               # Then I need to Drop the previous variables that will be replaced with new ones
                                               drop.cols <- c('fMi','beta0','beta1','p')
                                               assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",i))) %>% dplyr::select(-drop.cols))
                                               # Then Add the new betas and fMis to the dataframe 
                                               assign(paste0("my_data",j), eval(as.symbol(paste0("my_data",j))) %>% 
                                                        dplyr::mutate(beta0 = eval(as.symbol(paste0("nextbetas_",j)))[1], 
                                                                      beta1 = eval(as.symbol(paste0("nextbetas_",j)))[2], 
                                                                      iteration = j))
                                               mergeCols <- c("coverage","Y")
                                               assign(paste0("my_data",j,sep=""),dplyr::full_join(eval(as.symbol(paste0("nextfMis_",j))), eval(as.symbol(paste0("my_data",j))), by = mergeCols))
                                               #Now I have an updated dataframe of fMis and betas 
                                               
                                               #Next I need to Calculate the new P's for this 
                                               
                                               assign(paste0("my_data",j, sep = ""), calculate_p(eval(as.symbol(paste0("my_data",j)))))
                                               
                                               mydata_list[j] <- list(eval(as.symbol(paste0("my_data",j))))
                                               # Then calculate my LL for that iteration 
                                               LL_list[j] <- e_step_LL(eval(as.symbol(paste0("my_data",j))))
                                               
                                               abs_tol <- abs(LL_list[[j]] - LL_list[[i]])
                                               rel_tol <- (abs(LL_list[[j]] - LL_list[[i]]))/abs(LL_list[[i]])
                                               
                                               betas_next <- mydata_list[[j]] %>% dplyr::select(beta0,beta1) %>% distinct() %>% as.numeric()
                                               betas_previous <- mydata_list[[i]] %>% dplyr::select(beta0,beta1) %>% distinct() %>% as.numeric()
                                               
                                               beta0_change <- abs(betas_next[[1]]-betas_previous[[1]])
                                               beta1_change <- abs(betas_next[[2]]-betas_previous[[2]])
                                               
                                               
                                               i = i + 1
                                               print(j)
                                               if (j == 1000){
                                                 rel_tol <- -1 
                                                 LL_list[[j]] <- "max iterations reached"
                                               }
                                             }
                                             length(mydata_list) -> em_iteration
                                             
                                             mydata_list[[em_iteration]] %>% 
                                               dplyr::select(beta0,beta1) %>% 
                                               dplyr::distinct() %>% 
                                               as.numeric() -> my_em_betas
                                             
                                             cbind(my_em_betas[[1]], my_em_betas[[2]],LL_list[[em_iteration]], em_iteration[[1]]) -> final_LL_list[[k]]
                                           }
                                           
                                           do.call(rbind, final_LL_list) %>% 
                                             dplyr::as_tibble() %>% 
                                             dplyr::filter(!is.na(V1)) %>% 
                                             dplyr::distinct() -> results_tibble
                                           results_tibble[which.max(results_tibble$V3),] -> my_results 
                                           if (dim(my_results)[[1]] == 0){
                                             print("model didn't converge")
                                             my_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA)
                                           } else {
                                             print("model completed!")
                                           }
                                           
                                           output <- list("beta0" = my_results[[1]],
                                                          "beta1" = my_results[[2]], 
                                                          "LL" = my_results[[3]], 
                                                          "iterations" = my_results[[4]])
                                         })






proposed_parametric <- new_method("proposed_parametric_bootstrap", "proposed_parametric_bootstrap", 
                                  method = function(model,draw){
                                    restricted_function_noepsilon = function(par, my_data){
                                      e = 0 # e is parameter argument 1 
                                      vBeta = par[1]
                                      vGamma = par[2]
                                      Y0 <- my_data %>% dplyr::filter(Y == 0)
                                      Y1 <- my_data %>% dplyr::filter(Y == 1)
                                      
                                      Mi0 = as.matrix(Y0['coverage'])
                                      #Mi0 = log(Mi0)
                                      # add an intercept to the predictor variables
                                      #Mi0 = cbind(1, Mi0)
                                      mX0 = cbind(rep(1, nrow(Y0)))
                                      
                                      Mi1 = as.matrix(Y1['coverage'])
                                      #Mi1 = log(Mi1)
                                      #Mi1 = cbind(1, Mi1)
                                      # mX1 = as.matrix(Y1['host_bin'])
                                      mX1 = cbind(rep(1, nrow(Y1)))
                                      
                                      LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                      (1-(1-exp(-(Mi0 %*% vGamma))))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                        sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                  (1/(1+exp(-(mX1 %*% vBeta))))*(1-exp(-(Mi1 %*% vGamma)))))
                                      return(LL)
                                    }
                                    unrestricted_function_noepsilon = function(par, my_data){
                                      e = 0 # e is parameter argument 1 
                                      vBeta = c(par[1],par[2]) # my beta vector is parameter arguments 2 and 3 
                                      vGamma = par[3] # my gamma vector is parameter arguments 4
                                      Y0 <- my_data %>% dplyr::filter(Y == 0)
                                      Y1 <- my_data %>% dplyr::filter(Y == 1)
                                      
                                      Mi0 = as.matrix(Y0['coverage']) # I'm going to model the regular coverage
                                      # Mi0 = log(Mi0)
                                      # add an intercept to the predictor variables
                                      # Mi0 = cbind(1, Mi0)
                                      mX0 = as.matrix(Y0['X1'])
                                      mX0 = cbind(1, mX0)
                                      
                                      Mi1 = as.matrix(Y1['coverage'])
                                      # Mi1 = log(Mi1) # modeling the log coverage
                                      # Mi1 = cbind(1, Mi1)
                                      mX1 = as.matrix(Y1['X1'])
                                      mX1 = cbind(1, mX1)
                                      
                                      LL <- sum(log((1-e)*exp(-(mX0 %*% vBeta))/(1+exp(-(mX0 %*% vBeta))) + # when Y = 0 
                                                      (1-(1-exp(-(Mi0*vGamma))))*(1/(1+exp(-(mX0 %*% vBeta)))))) + 
                                        sum(log(e*exp(-(mX1 %*% vBeta))/(1+exp(-(mX1 %*% vBeta))) +  # when Y = 1
                                                  ((1-exp(-(Mi1*vGamma)))/(1+exp(-(mX1 %*% vBeta))))))
                                      return(LL)
                                    }
                                    
                                    
                                    Y0 <- draw %>% dplyr::filter(Y == 0)
                                    Y1 <- draw %>% dplyr::filter(Y == 1)
                                    Mi0 = as.matrix(Y0['coverage'])
                                    Mi0_mean = mean(Mi0)
                                    
                                    Mi1 = as.matrix(Y1['coverage'])
                                    Mi1_mean = mean(Mi1) 
                                    
                                    gamma_init <- (-1/Mi1_mean)*log((Mi0_mean)/(Mi0_mean+Mi1_mean))
                                    
                                    checkWarning <- base::tryCatch(mylogit <- stats::glm(Y ~ X1, data = draw, family = "binomial"), 
                                                                   warning = function(w){
                                                                     if(base::grepl("algorithm did not converge",as.character(w)) == "TRUE"){
                                                                       z <- 1}
                                                                   })
                                    
                                    if (length(checkWarning[[1]]) == 2){
                                      init_values <- c(mylogit$coefficients[[1]],mylogit$coefficients[[2]],gamma_init[[1]]) 
                                      init_values2 <- c(mylogit$coefficients[[1]],gamma_init[[1]])
                                      
                                      init_list <- list()
                                      init_list[[1]] <- init_values
                                      init_list[[2]] <- c(mylogit$coefficients[[1]]/3,mylogit$coefficients[[2]]/3,gamma_init[[1]]/3)
                                      init_list[[3]] <- c(mylogit$coefficients[[1]]*2,mylogit$coefficients[[2]]*2,gamma_init[[1]]*2)
                                      init_list[[4]] <- c(0.1,0.1,0.1)
                                      init_list[[5]] <- c(0.00001,0.00001,0.00001)
                                      init_list[[6]] <- c(0.0001,0.0001,0.0001)
                                      init_list[[7]] <- c(0,0,1)
                                      init_list[[8]] <- c(0,0,15)
                                      init_list[[9]] <- c(0,0,100)
                                      init_list[[10]] <- c(0,0,200)
                                      init_list[[11]] <- c(0,0,300)
                                      init_list[[12]] <- c(0,0,500)
                                      
                                      init_list2 <- list()
                                      init_list2[[1]] <- init_values2
                                      init_list2[[2]] <- c(mylogit$coefficients[[1]]/3,gamma_init[[1]]/3)
                                      init_list2[[3]] <- c(mylogit$coefficients[[1]]*2,gamma_init[[1]]*2)
                                      init_list2[[4]] <- c(0.1,0.1)
                                      init_list2[[5]] <- c(0.00001,0.00001)
                                      init_list2[[6]] <- c(0.0001,0.0001)
                                      init_list2[[7]] <- c(0,1)
                                      init_list2[[8]] <- c(0,15)
                                      init_list2[[9]] <- c(0,100)
                                      init_list2[[10]] <- c(0,200)
                                      init_list2[[11]] <- c(0,300)
                                      init_list2[[12]] <- c(0,500)
                                      
                                    } else {
                                      init_list <- list() 
                                      init_list[[1]] <- c(0.1,0.1,0.1)
                                      init_list[[2]] <- c(0.00001,0.00001,0.00001)
                                      init_list[[3]] <- c(0.0001,0.0001,0.0001)
                                      init_list[[4]] <- c(0.1,0.1,gamma_init[[1]])
                                      init_list[[5]] <- c(1,1,gamma_init[[1]]*2)
                                      init_list[[6]] <- c(0,0,1)
                                      init_list[[7]] <- c(0,0,15)
                                      init_list[[8]] <- c(0,0,100)
                                      init_list[[9]] <- c(0,0,200)
                                      init_list[[10]] <- c(0,0,300)
                                      init_list[[11]] <- c(0,0,500)
                                      
                                      init_list2 <- list()
                                      init_list2[[1]] <- c(0.1,0.1)
                                      init_list2[[2]] <- c(0.00001,0.00001)
                                      init_list2[[3]] <- c(0.0001,0.0001)
                                      init_list2[[4]] <- c(0.1,gamma_init[[1]])
                                      init_list2[[5]] <- c(1,gamma_init[[1]]*2)
                                      init_list2[[6]] <- c(0,1)
                                      init_list2[[7]] <- c(0,15)
                                      init_list2[[8]] <- c(0,100)
                                      init_list2[[9]] <- c(0,200)
                                      init_list2[[10]] <- c(0,300)
                                      init_list2[[11]] <- c(0,500)
                                    }
                                    
                                    restricted_list <- list()
                                    for (i in 1:length(init_list2)){  
                                      restricted_model = optimx::optimx(init_list2[[i]], restricted_function_noepsilon, my_data = draw, 
                                                                        method = 'L-BFGS-B',
                                                                        lower=c(-Inf,0),
                                                                        upper=c(Inf, Inf),
                                                                        control = list(fnscale= -1)) 
                                      cbind(restricted_model$p1, restricted_model$p2,restricted_model$value) -> restricted_list[[i]]
                                    }  
                                    do.call(rbind, restricted_list) %>% 
                                      dplyr::as_tibble() %>% 
                                      dplyr::filter(!is.na(V1)) %>% 
                                      dplyr::distinct() -> restricted_tibble
                                    restricted_tibble[which.max(restricted_tibble$V3),] -> restricted_results 
                                    if (dim(restricted_results)[[1]] == 0){
                                      print("restricted model no epsilon didn't converge")
                                      restricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA)
                                      LRT_hat = NA
                                      pvalue = NA
                                    } else {
                                      print("restricted model no epsilon completed!")
                                    }
                                    unrestricted_list <- list() 
                                    for (j in 1:length(init_list)){
                                      unrestricted_model = optimx::optimx(init_list[[j]], unrestricted_function_noepsilon, my_data = draw, 
                                                                          method = 'L-BFGS-B',
                                                                          lower= c(-Inf,-Inf, 0),
                                                                          upper=c(Inf, Inf, Inf),
                                                                          control = list(trace = 0, fnscale= -1))
                                      cbind(unrestricted_model$p1, unrestricted_model$p2,unrestricted_model$p3,unrestricted_model$value) -> unrestricted_list[[j]]
                                    }
                                    do.call(rbind, unrestricted_list) %>% 
                                      dplyr::as_tibble() %>% 
                                      dplyr::filter(!is.na(V1)) %>% 
                                      dplyr::distinct() -> unrestricted_tibble
                                    unrestricted_tibble[which.max(unrestricted_tibble$V4),] -> unrestricted_results 
                                    if (dim(unrestricted_results)[[1]] == 0){
                                      print("unrestricted model no epsilon didn't converge")
                                      unrestricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA)
                                      LRT_hat = NA 
                                      pvalue = NA
                                    } else {
                                      print("unrestricted model no epsilon completed!")
                                      LRT_hat = 2*abs(unrestricted_results$V4 - restricted_results$V3) }
                                    
                                    # Do I need to set a seed? I think the simulator should handle this... 
                                    if (!is.na(LRT_hat)){
                                      null_beta0 <- restricted_results$V1
                                      null_gamma1 <- restricted_results$V2 
                                      null_beta1 <- 0 
                                      X1 <- draw$X1
                                      coverage <- draw$coverage 
                                      LRT_list <- list()
                                      for (k in 1:1000){
                                        boot_beta.table <- cbind(X1,null_beta0,null_beta1) %>% dplyr::as_tibble()
                                        boot_beta.table$lambda_props = exp(boot_beta.table$null_beta0 + boot_beta.table$null_beta1*boot_beta.table$X1)/(1+exp(boot_beta.table$null_beta0+boot_beta.table$null_beta1*boot_beta.table$X1))
                                        # Draw from a Bernoulli distribution with p = lambda for both groups 
                                        boot_beta.table$lambda <- rbinom(nrow(draw), 1, boot_beta.table$lambda_props)
                                        # our gene's true presence is lambda
                                        boot_gamma.table <- cbind(boot_beta.table,coverage,null_gamma1) %>% dplyr::as_tibble()
                                        # gamma.table$lambda_bin <- ifelse(gamma.table$lambda >= 0.5, 1, 0)
                                        boot_gamma.table$Y_prop <- ifelse(boot_gamma.table$lambda == 1, 1-exp(-null_gamma1*boot_gamma.table$coverage), 0)
                                        boot_gamma.table$Y <- rbinom(nrow(draw), 1, boot_gamma.table$Y_prop)
                                        my_bootdata <- boot_gamma.table 
                                        
                                        Y0 <- my_bootdata %>% dplyr::filter(Y == 0)
                                        Y1 <- my_bootdata %>% dplyr::filter(Y == 1)
                                        Mi0 = as.matrix(Y0['coverage'])
                                        Mi0_mean = mean(Mi0)
                                        
                                        Mi1 = as.matrix(Y1['coverage'])
                                        Mi1_mean = mean(Mi1) 
                                        
                                        gamma_init <- (-1/Mi1_mean)*log((Mi0_mean)/(Mi0_mean+Mi1_mean))
                                        
                                        checkWarning <- base::tryCatch(mylogit <- stats::glm(Y ~ X1, data = my_bootdata, family = "binomial"), 
                                                                       warning = function(w){
                                                                         if(base::grepl("algorithm did not converge",as.character(w)) == "TRUE"){
                                                                           z <- 1}
                                                                       })
                                        
                                        if (length(checkWarning[[1]]) == 2){
                                          boot_init_values <- c(mylogit$coefficients[[1]],mylogit$coefficients[[2]],gamma_init[[1]]) 
                                          boot_init_values2 <- c(mylogit$coefficients[[1]],gamma_init[[1]])
                                          
                                          boot_init_list <- list()
                                          boot_init_list[[1]] <- boot_init_values
                                          boot_init_list[[2]] <- c(mylogit$coefficients[[1]]/3,mylogit$coefficients[[2]]/3,gamma_init[[1]]/3)
                                          boot_init_list[[3]] <- c(mylogit$coefficients[[1]]*2,mylogit$coefficients[[2]]*2,gamma_init[[1]]*2)
                                          boot_init_list[[4]] <- c(0.1,0.1,0.1)
                                          boot_init_list[[5]] <- c(0.00001,0.00001,0.00001)
                                          boot_init_list[[6]] <- c(0.0001,0.0001,0.0001)
                                          boot_init_list[[7]] <- c(0,0,1)
                                          boot_init_list[[8]] <- c(0,0,15)
                                          boot_init_list[[9]] <- c(0,0,100)
                                          boot_init_list[[10]] <- c(0,0,200)
                                          boot_init_list[[11]] <- c(0,0,300)
                                          boot_init_list[[12]] <- c(0,0,500)
                                          
                                          boot_init_list2 <- list()
                                          boot_init_list2[[1]] <- boot_init_values2
                                          boot_init_list2[[2]] <- c(mylogit$coefficients[[1]]/3,gamma_init[[1]]/3)
                                          boot_init_list2[[3]] <- c(mylogit$coefficients[[1]]*2,gamma_init[[1]]*2)
                                          boot_init_list2[[4]] <- c(0.1,0.1)
                                          boot_init_list2[[5]] <- c(0.00001,0.00001)
                                          boot_init_list2[[6]] <- c(0.0001,0.0001)
                                          boot_init_list2[[7]] <- c(0,1)
                                          boot_init_list2[[8]] <- c(0,15)
                                          boot_init_list2[[9]] <- c(0,100)
                                          boot_init_list2[[10]] <- c(0,200)
                                          boot_init_list2[[11]] <- c(0,300)
                                          boot_init_list2[[12]] <- c(0,500)
                                          
                                        } else {
                                          boot_init_list <- list() 
                                          boot_init_list[[1]] <- c(0.1,0.1,0.1)
                                          boot_init_list[[2]] <- c(0.00001,0.00001,0.00001)
                                          boot_init_list[[3]] <- c(0.0001,0.0001,0.0001)
                                          boot_init_list[[4]] <- c(0.1,0.1,gamma_init[[1]])
                                          boot_init_list[[5]] <- c(1,1,gamma_init[[1]]*2)
                                          boot_init_list[[6]] <- c(0,0,1)
                                          boot_init_list[[7]] <- c(0,0,15)
                                          boot_init_list[[8]] <- c(0,0,100)
                                          boot_init_list[[9]] <- c(0,0,200)
                                          boot_init_list[[10]] <- c(0,0,300)
                                          boot_init_list[[11]] <- c(0,0,500)
                                          
                                          boot_init_list2 <- list()
                                          boot_init_list2[[1]] <- c(0.1,0.1)
                                          boot_init_list2[[2]] <- c(0.00001,0.00001)
                                          boot_init_list2[[3]] <- c(0.0001,0.0001)
                                          boot_init_list2[[4]] <- c(0.1,gamma_init[[1]])
                                          boot_init_list2[[5]] <- c(1,gamma_init[[1]]*2)
                                          boot_init_list2[[6]] <- c(0,1)
                                          boot_init_list2[[7]] <- c(0,15)
                                          boot_init_list2[[8]] <- c(0,100)
                                          boot_init_list2[[9]] <- c(0,200)
                                          boot_init_list2[[10]] <- c(0,300)
                                          boot_init_list2[[11]] <- c(0,500)
                                        }
                                        
                                        boot_restricted_list <- list()
                                        for (l in 1:length(boot_init_list2)){  
                                          boot_restricted_model = optimx::optimx(boot_init_list2[[l]], restricted_function_noepsilon, my_data = my_bootdata, 
                                                                                 method = 'L-BFGS-B',
                                                                                 lower=c(-Inf,0),
                                                                                 upper=c(Inf, Inf),
                                                                                 control = list(fnscale= -1)) 
                                          cbind(boot_restricted_model$p1, boot_restricted_model$p2, boot_restricted_model$value) -> boot_restricted_list[[l]]
                                        }  
                                        do.call(rbind, boot_restricted_list) %>% 
                                          dplyr::as_tibble() %>% 
                                          dplyr::filter(!is.na(V1)) %>% 
                                          dplyr::distinct() -> boot_restricted_tibble
                                        boot_restricted_tibble[which.max(boot_restricted_tibble$V3),] -> boot_restricted_results 
                                        if (dim(boot_restricted_results)[[1]] == 0){
                                          print(paste("restricted model bootstrap",j, "didn't converge"))
                                          boot_restricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA)
                                          LRT_boot = NA
                                        } else {
                                          print(paste("restricted model bootstrap",j,"completed!"))
                                        }
                                        boot_unrestricted_list <- list() 
                                        for (m in 1:length(boot_init_list)){
                                          boot_unrestricted_model = optimx::optimx(boot_init_list[[m]], unrestricted_function_noepsilon, my_data = my_bootdata, 
                                                                                   method = 'L-BFGS-B',
                                                                                   lower= c(-Inf,-Inf, 0),
                                                                                   upper=c(Inf, Inf, Inf),
                                                                                   control = list(trace = 0, fnscale= -1))
                                          cbind(boot_unrestricted_model$p1, boot_unrestricted_model$p2, boot_unrestricted_model$p3, boot_unrestricted_model$value) -> boot_unrestricted_list[[m]]
                                        }
                                        do.call(rbind, boot_unrestricted_list) %>% 
                                          dplyr::as_tibble() %>% 
                                          dplyr::filter(!is.na(V1)) %>% 
                                          dplyr::distinct() -> boot_unrestricted_tibble
                                        boot_unrestricted_tibble[which.max(boot_unrestricted_tibble$V4),] -> boot_unrestricted_results 
                                        if (dim(boot_unrestricted_results)[[1]] == 0){
                                          print(paste("unrestricted model bootstrap",k,"didn't converge"))
                                          boot_unrestricted_results <- dplyr::tibble(V1 = NA, V2 = NA, V3 = NA, V4 = NA)
                                          LRT_boot = NA 
                                          pvalue = NA
                                        } else {
                                          print(paste("unrestricted model bootstrap",k,"completed!"))
                                          LRT_boot <- 2*abs(boot_unrestricted_results$V4 - boot_restricted_results$V3)}
                                        LRT_list[[k]] <- LRT_boot
                                      }
                                      do.call(rbind, LRT_list) %>% 
                                        dplyr::as_tibble() %>% 
                                        dplyr::filter(!is.na(V1)) %>% 
                                        cbind(LRT_hat) %>% 
                                        dplyr::mutate(Q = ifelse(V1 >= LRT_hat,1,0)) -> LRT_table 
                                      pvalue = (1/(nrow(LRT_table)+1))*(1 + sum(LRT_table$Q))
                                    }else {pvalue = NA}
                                    output <- list("pvalue" = pvalue)
                                  }) 



GLM_Rao <- new_method("GLM_Rao","GLM_Rao", 
                      method = function(model,draw){
                        if(length(draw$Y %>% unique()) == 1){
                          pvalue <- 1
                          Rao <- NA
                          beta0 <- NA 
                          beta1 <- NA
                        } else if(length(draw$X1 %>% unique()) == 1) {
                          pvalue <- 1
                          Rao <- NA
                          beta0 <- NA 
                          beta1 <- NA
                          
                        } else {
                          draw %>% 
                            stats::glm(Y ~ X1,.,
                                       family = binomial(link = "logit")) -> naive_glm
                          naive_glm %>%
                            anova(test="Rao") -> anova_output 
                          
                          anova_output$Rao[2] -> Rao
                          anova_output$`Pr(>Chi)`[2] -> pvalue 
                          naive_glm$coefficients[[1]] -> beta0
                          naive_glm$coefficients[[2]] -> beta1 
                        }
                        
                        output <- list("Rao" = Rao,
                                       "pvalue" = pvalue,
                                       "beta0" = beta0,
                                       "beta1" = beta1)
                      })

GLM_LRT <- new_method("GLM_LRT","GLM_LRT", 
                      method = function(model,draw){
                        if(length(draw$Y %>% unique()) == 1){
                          pvalue <- 1
                          beta0 <- NA 
                          beta1 <- NA
                        } else if(length(draw$X1 %>% unique()) == 1) {
                          pvalue <- 1
                          beta0 <- NA 
                          beta1 <- NA
                        } else {
                          draw %>% 
                            stats::glm(Y ~ X1,.,
                                       family = binomial(link = "logit")) -> full_model
                          draw %>% 
                            stats::glm(Y ~ 1,., family = binomial(link="logit")) -> reduced_model
                          
                          anova(full_model, reduced_model, test = "Chisq") -> anova_output 
                          
                          anova_output$`Pr(>Chi)`[2] -> pvalue 
                          full_model$coefficients[[1]] -> beta0
                          full_model$coefficients[[2]] -> beta1
                        }
                        
                        output <- list("pvalue" = pvalue,
                                       "beta0" = beta0,
                                       "beta1" = beta1)
                      })


Fisher <- new_method("FisherTest","FisherTest", 
                     method = function(model,draw){
                       if(length(draw$Y %>% unique()) == 1){
                         pvalue <- 1
                         beta0 <- NA 
                         beta1 <- NA
                       } else {
                         table(draw$Y,draw$X1) -> contingencytable
                         fisher.test(contingencytable) -> fisher_results
                         fisher_results$p.value -> pvalue
                       }
                       output <- list("pvalue" = pvalue)
                     })