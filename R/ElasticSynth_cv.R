#' ElasticSynth
#'
#' A synthetic controls method for comparitive case studies
#' that utilizes elastic net regression -- glmnet -- to assign unit weights. 
#' Utilizes a cross-validation mechanism for selecting optimal L1 and L2 regularization parameters, as well as lambda.
#' 
#' @param OutcomeMatrix a t x n matrix where each row is the outcome variable of interest at time t for each unit
#' @param treated the column index for the treated unit in both PredictorMatrix and OutcomeMatrix. For multiple treated units, see 'ElasticSynthRun'
#' @param pre a vector of row indices indicating the pre period (Ex. 1:40). This is the period over which the algorithm will optimize weights.
#' @param post a vector of row indices indicating the post period (Ex. 41:50)
#' @param alphas glmnet param - a vector of values for alpha to be fitted in cross validation
#' @param lambdas glmnet param - a vector of values for lambda to be fitted in cross validation
#' @param start_cv which time period to start Hyndman's cv method
#' @param end_cv which time period to end Hyndman's cv method (usually last pre period)
#' @param cv_step number of time periods to forecast out
#' @param lower_limit_weights The lower limit value of weight that can be placed on any unit. Default is zero (non-negative constraint). Change if desired, but be wary of overfitting.
#' @param upper_limit_weights The upper limit value of weight that can be placed on any unit. Default is one.
#' @param placebo If TRUE, run cv for all units. If false, run cv for only treated unit
#' @param verbose Print unit status during cross validation?
#' @return optimal alpha and lambda for each unit as well as pre-period cv mse
#' @export

ElasticSynth_cv = function(
  OutcomeMatrix, 
  treated = 1, 
  pre, 
  post,   
  alphas,
  lambdas,
  start_cv,
  end_cv,
  cv_step,
  lower_limit_weights = 0,
  upper_limit_weights = Inf,
  placebo = T,
  verbose = T) 
{
  
  
  
  suppressMessages(library(glmnet))
  suppressMessages(library(ggplot2))
  suppressMessages(library(lubridate))
  suppressMessages(library(dplyr))
  suppressMessages(library(tidyr))
  suppressMessages(library(Metrics))
  suppressMessages(library(reshape2))
  
  set.seed(12345)
  
  ####  data structure
  
  # Time Periods for units minus treated
  Y0 = as.matrix(OutcomeMatrix[,-treated])
  # Time Periods for treated
  Y1 = as.matrix(OutcomeMatrix[,treated, drop = F])
  # Combined
  Y  = as.matrix(cbind(Y1, Y0))
  # Pre period units minus treated
  Z0 = as.matrix(Y0[pre,])
  # Pre period for treated
  Z1 = as.matrix(Y1[pre,1,drop = F])
  # Combined
  Z  = as.matrix(cbind(Z1, Z0))
  
  ##### Params
  
  # Number of units
  if (placebo) {
    N  = dim(Y)[2]
  }
  
  if (!placebo) {
    N = 1
  }
  # Total time periods
  Time = length(c(pre,post))
  # Pre period
  T0   = length(pre)
  # Post Period
  T1   = length(post)
  
  ## Hyndman CV fold structure
  number_folds =  end_cv - (start_cv - 1) - (step - 1)
  
  folds_ix     = list()
  for (i in 1:number_folds) {
    folds_ix[[i]] = (start_cv + (i - 1)):(start_cv + (i - 1) + (cv_step - 1))
  }
  
  
  ## Month frame for plotting
  na  = length(alphas)
  results = data.frame(unit = colnames(Y), alpha = rep(0, dim(Y)[2]), lambda = rep(0, dim(Y)[2]), cv_mse = rep(0, dim(Y)[2]))
  
  cat('\n\n\n\n******* Cross Validation *******n\n\n\n')
  
  for (i in 1:N) {
    
    if (verbose == T) {
      cat('*** Unit', colnames(Y[,i,drop = F]), '***\n\n')
    }
    
    Y1       = as.matrix(Y[,i, drop = F])
    unitName = gsub('_[0-9]+', '', colnames(Y1))
    Y0       = as.matrix(Y[,-c(1, i, grep(paste(unitName, '_', sep = ''), colnames(Y)))])
    Z1       = as.matrix(Z[,i, drop = F])
    Z0       = as.matrix(Z[,-c(1, i, grep(paste(unitName, '_', sep = ''), colnames(Z)))])
    cv_grid = data.frame()
    
    for (r in 1:number_folds) {
      # hyndman cv structure
      Z0_train = Z0[1:(min(folds_ix[[r]]) - 1),] 
      Z1_train = Z1[1:(min(folds_ix[[r]]) - 1),]
      Z0_test  = Z0[min(folds_ix[[r]]):max(folds_ix[[r]]),]
      Z1_test  = Z1[min(folds_ix[[r]]):max(folds_ix[[r]]),]
      alpha_grid = data.frame()
      for (j in 1:na) {
        a        = alphas[j]
        fit      = glmnet(x = Z0_train, 
                          y = Z1_train,
                          family = 'gaussian',
                          alpha = a,
                          lambda = lambdas,
                          standardize = FALSE,
                          intercept = FALSE,
                          lower.limits = lower_limit_weights,
                          upper.limits = upper_limit_weights)
        y_pred  = predict(fit, newx = Z0_test, type = 'response', s = fit$lambda)
        y_true  = Z1_test
        mse_lambda = as.numeric(apply(y_pred, 2, function(y_pred) Metrics::mse(y_true, y_pred)))
        lambda_ix  = which(mse_lambda == min(mse_lambda))
        alpha_grid_tmp = data.frame(fold = rep(r, length(lambdas)), alpha = rep(a, length(lambdas)), lambda = lambdas, mse = mse_lambda)
        alpha_grid = rbind(alpha_grid, alpha_grid_tmp)
      }
      cv_grid = rbind(cv_grid, alpha_grid)
    }
    final_grid = dcast(cv_grid, lambda + alpha ~ fold, value.var = 'mse')
    row_means  = apply(final_grid[,3:(number_folds+2)],1,mean)
    best_ix    = which(row_means == min(row_means))
    results$alpha[i] = final_grid[best_ix,'alpha']
    results$lambda[i] = final_grid[best_ix,'lambda']
    results$cv_mse[i] = row_means[best_ix]
  }
  
  
  # treated unit
  
  
  return(results)
  
  
}
