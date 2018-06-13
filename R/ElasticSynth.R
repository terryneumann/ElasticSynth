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
#' @param a_grid glmnet param - a vector of values for alpha to be fitted in cross validation
#' @param start_month the start month of the case study (first pre-period). This is used for plotting
#' @param end_month the end month of the case study (last post-period). This is used for plotting
#' @param time_unit The gap between observations in the case study. No default. A character string, containing one of "day", "week", "month", "quarter" or "year". This can optionally be preceded by an integer and a space, or followed by "s". Ex '3 months' or 'quarter'
#' @param penaltyMat glmnet param - (optional) a n-1 x n matrix of additional penalties placed upon glmnet for selecting certain units. Useful in the case of multiple treated units. Each column of the matrix represents the penalties to be placed on the control units. The treated unit is excluded from the column. See glmnet param penalty.factor for more details. Default is no additional penalties on any unit.
#' @param max_number_units glmnet param - The maximum number of units glmnet will select. Useful for forcing parsimonious models. Default allows all units to receive weight.
#' @param lower_limit_weights The lower limit value of weight that can be placed on any unit. Default is zero (non-negative constraint). Change if desired, but be wary of overfitting.
#' @param upper_limit_weights The upper limit value of weight that can be placed on any unit. Default is one.
#' @param nfolds The number of folds used in cross validation. The default is 5.
#' @param test Hypothesis test - 'lower', 'upper', or 'two-sided'
#' @param verbose Print unit status during cross validation?
#' @return list containing output weights for treated unit (w_final), the actual outcome values (Y_true), the fitted outcome values (Y_elast), optimal value of lambda (lambda_opt), optimal value of alpha, (alpha_opt), a dataframe of the results of placebo test (placebo_frame), a plot of the path of the treated vs. actual unit (path.plot), the deviation ratio of the fitted series (dev.ratio), and a plot of the placebo test results for the treated unit. 
#' @export

ElasticSynth = function(
  OutcomeMatrix, 
  treated, 
  pre, 
  post,   
  a_grid =seq(from = 0.1, to = 0.9, by = 0.1),
  start_month,
  end_month,
  time_unit,
  max_number_units = ncol(OutcomeMatrix),
  lower_limit_weights = 0,
  upper_limit_weights = 1,
  nfolds = 5,
  test,
  err_alpha_lambda_opt = NULL,
  verbose = T) 
{
  
  
  suppressMessages(library(glmnet))
  suppressMessages(library(ggplot2))
  suppressMessages(library(lubridate))
  suppressMessages(library(dplyr))
  suppressMessages(library(tidyr))
  suppressMessages(library(Metrics))
  
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
  N  = dim(Y)[2]
  # Total time periods
  Time = length(c(pre,post))
  # Pre period
  T0   = length(pre)
  # Post Period
  T1   = length(post)
  
  
  ## Month frame for plotting
  start_month = as.Date(start_month)
  end_month = as.Date(end_month)
  month_seq  = seq(start_month, end_month, by = time_unit)
  month_join = data.frame(time = 1:max(post), month = month_seq)
  na  = length(a_grid)
  
  
  cat('\n\n\n\n******* Cross Validation *******n\n\n\n')
  
  if (is.null(err_alpha_lambda_opt)) {
    err_alpha_lambda    = expand.grid(a = a_grid, opt_lambda = 0, error = 0, unit = colnames(Y))
    
    
    for (i in 1:N) {
      
      if (verbose == T) {
        cat('*** Unit', colnames(Y[,i,drop = F]), '***\n\n')
      }
      
      for (j in 1:na) {
        a       = a_grid[j]
        Y1       = as.matrix(Y[,i, drop = F])
        unitName = gsub('_[0-9]+', '', colnames(Y1))
        Y0       = as.matrix(Y[,-c(1, i, grep(paste(unitName, '_', sep = ''), colnames(Y)))])
        Z1       = as.matrix(Z[,i, drop = F])
        Z0       = as.matrix(Z[,-c(1, i, grep(paste(unitName, '_', sep = ''), colnames(Z)))])
        V1      = scale(Z1, scale = FALSE)
        V0      = scale(Z0, scale = FALSE)
        fit     = cv.glmnet(x = V0, y = V1,
                            alpha = a,
                            standardize = FALSE,
                            intercept = FALSE,
                            lower.limits = lower_limit_weights,
                            upper.limits = upper_limit_weights,
                            pmax = max_number_units,
                            nfolds = nfolds,
                            type.measure = 'mse')
        err_alpha_lambda$error[j + (i-1)*na]      = fit$cvm[fit$lambda == fit$lambda.min]
        err_alpha_lambda$opt_lambda[j + (i-1)*na] = fit$lambda[fit$lambda == fit$lambda.min]
        
      }
    }
    
    err_alpha_lambda_opt = err_alpha_lambda %>% 
      group_by(unit) %>%
      summarise(a = a[which.min(error)], lambda = opt_lambda[which.min(error)], error = min(error))
    
  }
  else {
    err_alpha_lambda    = expand.grid(a = a_grid, opt_lambda = 0, error = 0, unit = colnames(Y)[1])
    for (j in 1:na) {
      a        = a_grid[j]
      Y1       = as.matrix(Y[,1, drop = F])
      unitName = gsub('_[0-9]+', '', colnames(Y1))
      Y0       = as.matrix(Y[,-c(1)])
      Z1       = as.matrix(Z[,i, drop = F])
      Z0       = as.matrix(Z[,-c(1)])
      V1      = scale(Z1, scale = FALSE)
      V0      = scale(Z0, scale = FALSE)
      fit     = cv.glmnet(x = V0, y = V1,
                          alpha = a,
                          standardize = FALSE,
                          intercept = FALSE,
                          lower.limits = lower_limit_weights,
                          upper.limits = upper_limit_weights,
                          pmax = max_number_units,
                          nfolds = nfolds)
      err_alpha_lambda$error[j]      = fit$cvm[fit$lambda == fit$lambda.min]
      err_alpha_lambda$opt_lambda[j] = fit$lambda[fit$lambda == fit$lambda.min]
      
    }
    
    err_alpha_lambda_opt_new = err_alpha_lambda %>% 
      group_by(unit) %>%
      summarise(a = a[which.min(error)], lambda = opt_lambda[which.min(error)], error = min(error))
    
    err_alpha_lambda_opt = err_alpha_lambda_opt[-1,]
    err_alpha_lambda_opt = rbind(err_alpha_lambda_opt_new, err_alpha_lambda_opt)
  }
  ###### Cross validation for many lambda over given alpha
  
  err_alpha_lambda_opt$error_post = rep(0, N)
  int_elast   = matrix(0, nrow = 1, ncol = 1)
  Y_elast     = matrix(0, nrow = Time, ncol = 1)
  Y_true      = matrix(0, nrow = Time, ncol = 1)
  old_frame   = data.frame()
  for (i in 1:N) {
    if (i != 1) {
      old_frame = new_frame
    }
    Y1 = as.matrix(Y[,i,drop = F])
    unitName = gsub('_[0-9]+', '', colnames(Y1))
    Y0 = as.matrix(Y[,-c(i, grep(paste(unitName, '_', sep = ''), colnames(Y)))])
    Z1 = as.matrix(Z[,i, drop = F])
    Z0 = as.matrix(Z[,-c(i, grep(paste(unitName, '_', sep = ''), colnames(Y)))])
    V1 = rbind(scale(Z1, scale = FALSE))
    V0 = rbind(scale(Z0, scale = FALSE))
    
    # Fit elast
    fit_final   = glmnet(x = V0, y = V1,
                         alpha = err_alpha_lambda_opt$a[i],
                         lambda = err_alpha_lambda_opt$lambda[i],
                         standardize = FALSE,
                         intercept = FALSE,
                         pmax = max_number_units,
                         lower.limits = lower_limit_weights,
                         upper.limits = upper_limit_weights)
    w           = as.matrix(coef(fit_final, s = err_alpha_lambda_opt$lambda[i]))[-1,]
    if (i == 1) {

      w_final     = w
      int_elast   = as.matrix(apply(Z1 - Z0 %*% w_final, 2, mean))
      Y_elast     = int_elast[rep(1, Time),] + Y0 %*% w_final
      Y_elast_final = Y_elast
      Y_true      = Y1[c(1:max(pre),post)]
      Y_true_final = Y_true
      gaps          = Y_true[c(1:max(pre),post)] - Y_elast[c(1:max(pre),post)]
      Y_smape       = Metrics::smape(Y_elast[1:max(pre)], Y_true[1:max(pre)])
      
      
      plotFrame   = data.frame(time = c(pre,post), Y_true = Y_true, Y_elast = Y_elast, gaps = gaps)
      plotFrame   = merge(plotFrame, month_join, by = 'time', all.x = TRUE)
      min_post_month = plotFrame$month[plotFrame$time == min(post)]
      max_post_month = plotFrame$month[plotFrame$time == max(post)]
      # 
      path.plot   = ggplot(plotFrame) +
        geom_line(aes(x = as.Date(month), y = Y_true), colour = 'black', size = 2) +
        geom_line(aes(x = as.Date(month), y = Y_elast), colour = 'indianred3', linetype = 'longdash', size = 2) +
        geom_rect(aes(xmin = min_post_month, xmax = max_post_month, ymin = -Inf, ymax = Inf), color="olivedrab2",
                  alpha=0.01,
                  inherit.aes = FALSE) +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"),
              axis.line.x = element_line(colour = 'black'),
              axis.line.y = element_line(colour = 'black')) +
        xlab('Time Period') +
        ylab('Outcome') +  
        xlim(c(as.Date(start_month), as.Date(end_month %m+% months(2))))
      
      
      
    }
    
    int_elast   = as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
    Y_elast     = int_elast[rep(1, Time),] + Y0 %*% w
    Y_true      = Y1
    err_alpha_lambda_opt$error_post[i] = Metrics::mse(Y_true[post], Y_elast[post])
    gaps        = Y_true[c(1:max(pre),post)] - Y_elast[c(1:max(pre),post)]
    new_frame   = data.frame(gaps = unlist(gaps),
                             y_true = Y_true[c(1:max(pre),post)], 
                             y_elast = Y_elast[c(1:max(pre),post)],
                             time = c(1:max(pre),post), 
                             unit_type = ifelse(i == 1, 'Treated Unit', 'Control Unit Distribution'), 
                             unit = colnames(Y)[i])
    new_frame   = rbind(old_frame, new_frame)
    
  }
  new_frame     = merge(new_frame, month_join, by = 'time', all.x = TRUE)
  percent_frame = subset(new_frame, time %in% post) %>%
    mutate(p_stat_lower = rank(gaps)/length(gaps),
           p_stat_above = 1 - (rank(gaps)/length(gaps)),
           p_stat_two_tail = ifelse(rank(gaps)/length(gaps) > 0.5, 1 - (rank(gaps)/length(gaps)), rank(gaps)/length(gaps)))
  
  if (test == 'lower') {
    p_stat = subset(percent_frame, unit_type == 'Treated Unit')$p_stat_lower
  }
  else if (test == 'upper') {
    p_stat = subset(percent_frame, unit_type == 'Treated Unit')$p_stat_above
  }
  else if (test == 'two-tail') {
    p_stat = subset(percent_frame, unit_type == 'Treated Unit')$p_stat_two_tail
  }
  
  
  min_post_month = unique(new_frame$month[new_frame$time == min(post)])
  max_post_month = unique(new_frame$month[new_frame$time == max(post)])
  far_x_axis_month = max_post_month %m+% months(2)
  near_x_axis_month = max_post_month %m-% months(3*10)
  
  
  placebo = ggplot() +
    geom_line(data = subset(new_frame, unit_type == 'Control Unit Distribution'), aes(x = month, y = gaps, group = unit),colour = 'grey', size = 0.5, alpha = 0.4) +
    geom_line(data = subset(new_frame, unit_type == 'Treated Unit'), aes(x = month, y = gaps), size = 2, colour = 'black') +
    geom_hline(aes(yintercept = 0)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"), 
          axis.line.x = element_line(colour = 'black'),
          axis.line.y = element_line(colour = 'black')) +
    geom_rect(aes(xmin = min_post_month, xmax = max_post_month, ymin = -Inf, ymax = Inf), color="black",
              alpha=0.01,
              inherit.aes = FALSE) +
    ggtitle('Gaps Between Actual Unit and Synthetic Unit,\nwith Results from Placebo Test in Post Period') + 
    xlim(c(as.Date(near_x_axis_month), as.Date(far_x_axis_month)))
  
  
  
  # treated unit
  
  
  return(list(w = w_final, 
              Y_true = Y_true_final, 
              Y_elast = Y_elast_final, 
              fit = fit_final, 
              #    lambda_opt = lambda_opt, 
              #    alpha_opt = a_opt,
              placebo_frame = new_frame,
              placebo = placebo,
              path.plot = path.plot,
              smape = Y_smape,
              dev.ratio = dev.ratio,
              err_alpha_lambda_opt = err_alpha_lambda_opt,
              p_stat = p_stat,
              mse_frame = err_alpha_lambda_opt))
  
  
}
