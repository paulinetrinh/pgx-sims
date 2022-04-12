## @knitr evals
pvalue <- new_metric("pvalue", "pvalue",
                     metric = function(model,out){
                       return(out$pvalue)
                     })


gamma1_estimate <- new_metric("gamma1_hat", "gamma1_hat", 
                              metric = function(model,out){
                                return(out$gamma1)
                              })

beta0_estimate <- new_metric("beta0_hat", "beta0_hat", 
                             metric = function(model,out){
                               return(out$beta0)
                             })

beta1_estimate <- new_metric("beta1_hat", "beta1_hat", 
                             metric = function(model,out){
                               return(out$beta1)
                             })

iteration_eval <- new_metric("iteration", "iteration", 
                             metric = function(model,out){
                               return(out$iteration)
                             })

ftilde_eval <- new_metric("ftilde_change", "ftilde_change", 
                             metric = function(model,out){
                               return(out$ftilde_change)
                             })

MSE_beta1 <- new_metric("MSE_beta1", "MSE_beta1", 
                        metric = function(model,out){
                          
                        }
                        )

LL_unrestricted <- new_metric("LL_unrestricted", "LL_unrestricted", 
                             metric = function(model,out){
                               return(out$LL_unrestricted)
                             })

LL_restricted <- new_metric("LL_restricted", "LL_restricted", 
                            metric = function(model,out){
                              return(out$LL_restricted)
                            })
pvalue_noabs <- new_metric("pvalue_noabs", "pvalue_noabs",
                           metric = function(model,out){
                             pvalue_noabs <- stats::pchisq(out$LRT_noabs, df = 1, lower.tail = FALSE)
                             return(pvalue_noabs)
                           })
LRT_noabs <- new_metric("LRT_noabs", "LRT_noabs", 
                       metric = function(model,out){
                         return(out$LRT_noabs)
                       })

LRT_abs <- new_metric("LRT_abs", "LRT_abs", 
                       metric = function(model,out){
                         return(out$LRT_abs)
                       })