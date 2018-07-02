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
#' @param test Hypothesis test - 'lower', 'upper', or 'two-sided'
#' @param verbose Print unit status during cross validation?
#' @return list containing output weights for treated unit (w_final), the actual outcome values (Y_true), the fitted outcome values (Y_elast), optimal value of lambda (lambda_opt), optimal value of alpha, (alpha_opt), a dataframe of the results of placebo test (placebo_frame), a plot of the path of the treated vs. actual unit (path.plot), the deviation ratio of the fitted series (dev.ratio), and a plot of the placebo test results for the treated unit. 
#' @export

ElasticSynth_fit_subunit = function(
  treated_units,
  sub_units,
  placebo_units,
  treated, 
  pre, 
  post,
  test,
  start_month,
  end_month,
  time_unit,
  lower_limit_weights = 0,
  upper_limit_weights = 1,
  cv_results = NULL,
  storage_directory = NULL,
  verbose = T,
  placebo_graphs = T) 

{
  
  
  
  ####  data structure
  OutcomeMatrix = cbind(treated_units, placebo_units)
  # Time Periods for units minus treated
  Y0 = as.matrix(OutcomeMatrix[,-treated])
  # Time Periods for treated
  Y1 = as.matrix(OutcomeMatrix[,treated, drop = F])
  # Combined
  Y  = as.matrix(cbind(Y1, Y0))
  # Pre period units minus treated
  # Pre period for treated
  Z1 = as.matrix(Y1[pre,1,drop = F])
  # Combined
  Z  = as.matrix(cbind(Z1, Y0[pre,]))
  
  ##### Params
  n_tr = dim(treated_units)[2]
  # Number of units

  # Total time periods
  Time = length(c(pre,post))
  # Pre period
  T0   = length(pre)
  # Post Period
  T1   = length(post)
  
  N = dim(OutcomeMatrix)[2]
  
  
  cv_results$error_post = rep(0, N)
  old_frame   = data.frame()
  start_month = as.Date(start_month)
  end_month = as.Date(end_month)
  month_seq  = seq(start_month, end_month, by = time_unit)
  month_join = data.frame(time = min(pre):max(post), month = month_seq)
  w_final = list(list())
  
  for (i in 1:N) {
    

    if (i != 1) {
      old_frame = new_frame
    }
    
    Y1 = as.matrix(Y[,i,drop = F])
    unitName = sub('District ', '', names(OutcomeMatrix)[i])
    beats_minus_family = sub_units[c(min(pre):max(pre), post),grepl(substr(unitName, 1, 2), substr(names(sub_units), 6, 7)) == F,drop=F]
    Z1 = as.matrix(Z[,i, drop = F])
    Z0 = as.matrix(beats_minus_family[1:length(pre),])

    # Fit elast
    fit_final   = glmnet(x = Z0, 
                         y = Z1,
                         family = 'gaussian',
                         alpha = cv_results$alpha[i],
                         lambda = cv_results$lambda[i],
                         standardize = FALSE,
                         intercept = FALSE,
                         lower.limits = lower_limit_weights,
                         upper.limits = upper_limit_weights)
    
    w           = as.matrix(coef(fit_final, s = cv_results$lambda[i]))[-1,]
    
    if (i %in% 1:n_tr) {

      w_final[[i]]     = w
      int_elast        = as.matrix(apply(Z1 - Z0 %*% w_final[[i]], 2, mean))
      Y_elast          = int_elast[rep(1, Time),] + as.matrix(beats_minus_family) %*% w_final[[i]]
      Y_elast_final    = Y_elast
      Y_true           = Y1[c(min(pre):max(pre),post)]
      Y_true_final     = Y_true
      gaps             = Y_true[c(min(pre):max(pre),post)] - Y_elast[c(min(pre):max(pre),post)]
      
      
      
      plotFrame   = data.frame(time = c(min(pre):max(pre),post), Y_true = Y_true, Y_elast = Y_elast, gaps = gaps)
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
        xlim(c(as.Date(start_month), as.Date(end_month %m+% months(2)))) +
        ggtitle(paste('Actual (Black) vs Synthetic for', names(treated_units)[i], sep = ' '))
      
        ggsave(path.plot, filename = paste(storage_directory, 'path_', names(treated_units)[i], '.png', sep = ''), width = 12, height = 8)
      
      
    }
    
    int_elast                   = as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
    Y_elast                     = int_elast[rep(1, Time),] + as.matrix(beats_minus_family) %*% w
    Y_true                      = Y1
    gaps                        = Y_true[c(min(pre):max(pre),post)] - Y_elast[c(min(pre):max(pre),post)]
    cv_results$fitted_error_pre[i] = Metrics::mse(Y_true[pre], Y_elast[pre])
    cv_results$fitted_mape_pre[i]  = Metrics::mape(Y_true[pre], Y_elast[pre])
    cv_results$error_post[i]       = Metrics::mse(Y_true[post], Y_elast[post])
    cv_results$fitted_mape_post[i] = Metrics::mape(Y_true[post], Y_elast[post])
    cv_results$sign[i]             = ifelse(sum(gaps)<0, 'negative', 'postive')
    new_frame   = data.frame(gaps = unlist(gaps),
                             y_true = Y_true[c(min(pre):max(pre),post)], 
                             y_elast = Y_elast[c(min(pre):max(pre),post)],
                             time = c(min(pre):max(pre),post), 
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
  if (test == 'upper') {
    p_stat = subset(percent_frame, unit_type == 'Treated Unit')$p_stat_above
  }
  if (test == 'two-tail') {
    p_stat = subset(percent_frame, unit_type == 'Treated Unit')$p_stat_two_tail
  }
  
  
  min_post_month = unique(new_frame$month[new_frame$time == min(post)])
  max_post_month = unique(new_frame$month[new_frame$time == max(post)])
  far_x_axis_month = max_post_month %m+% months(2)
  near_x_axis_month = max_post_month %m-% months(3*10)
  
  if (placebo_graphs == T) {
    
    for (i in 1:n_tr) {
      
      placebo = ggplot() +
        geom_line(data = subset(new_frame, unit_type == 'Control Unit Distribution' & !(unit %in% names(treated_units))), aes(x = month, y = gaps, group = unit),colour = 'grey', size = 0.5, alpha = 0.4) +
        geom_line(data = subset(new_frame, unit == names(treated_units)[i]), aes(x = month, y = gaps), size = 2, colour = 'black') +
        geom_hline(aes(yintercept = 0)) +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"), 
              axis.line.x = element_line(colour = 'black'),
              axis.line.y = element_line(colour = 'black')) +
        geom_rect(aes(xmin = min_post_month, xmax = max_post_month, ymin = -Inf, ymax = Inf), color="black",
                  alpha=0.01,
                  inherit.aes = FALSE) +
        ggtitle(paste('Gaps Between Actual Unit and Synthetic Unit for ', names(treated_units)[i],'\nwith Results from Placebo Test in Post Period', sep = '')) + 
        xlim(c(as.Date(near_x_axis_month), as.Date(far_x_axis_month)))
      
      ggsave(placebo, filename = paste(storage_directory, 'placebo_', names(treated_units)[i], '.png', sep = ''), width = 12, height = 8)
      
      
    }
    
  }
  
  
  cv_results$cv_mse_ratio = cv_results$error_post/cv_results$cv_mse
  cv_results$cv_mse_ratio_p_value = 1 - rank(cv_results$cv_mse_ratio)/nrow(cv_results)
  cv_results$fitted_mse_ratio = cv_results$error_post/cv_results$fitted_error_pre
  cv_results$fitted_mape_ratio = cv_results$fitted_mape_post/cv_results$fitted_mape_pre
  
  ratio.plot = ggplot() +
    geom_point(data = subset(cv_results, unit %in% names(treated_units)), aes(x = cv_mse_ratio, y = 0.25, colour = sign), size = 3) +
    geom_segment(data = subset(cv_results, unit %in% names(treated_units)), aes(x = cv_mse_ratio, y = 0.25, yend = 0, xend = cv_mse_ratio, colour = sign), size = 1.1) +
    geom_density(data = subset(cv_results, !(unit %in% names(treated_units)) & cv_mse_ratio < quantile(cv_results$cv_mse_ratio, .995)), aes(x = cv_mse_ratio, y = ..density..), fill = 'grey', alpha = 0.6) +
    geom_text(data = subset(cv_results, unit %in% names(treated_units)), aes(x = cv_mse_ratio, y =.25, label = paste(unit, ': p = ', round(cv_mse_ratio_p_value, 3))), angle = 270, vjust = -.5) +
    scale_color_manual('Treatment Effect', values = c('blue','red')) +
    ggthemes::theme_gdocs() +
    scale_x_continuous(breaks = 1:round(quantile(cv_results$cv_mse_ratio, .995))) +
    ylab('Density') +
    ggtitle('Distribution of Test Statistic (Post Period MSE / Pre-period CV MSE) for control units\nwith Value and Effect for Treated Units Marked') +
    xlab('Post Period MSE / Pre Period CV MSE')
  
  ggsave(ratio.plot, filename = paste(storage_directory, 'ratio_plot.png', sep = ''), width = 12, height = 8)
  
  pre_period_cv_mse = ggplot(subset(cv_results, unit %in% names(OutcomeMatrix)[1:n_tr])) +
                        geom_bar(stat = 'identity', aes(x = unit, y = cv_mse, fill = unit), colour = 'black') +
                        coord_flip() +
                        ggthemes::theme_few() +
                        ggtitle('Pre-Period (Hyndman) Cross Validation MSE for Treated Units') +
                        ylab('Cross-Validation MSE') +
                        theme(text = element_text(size = 16),legend.position = 'none' )
  ggsave(pre_period_cv_mse, filename = paste(storage_directory, 'pre_period_cv_mse.png', sep = ''), width = 12, height = 8)
  
  
  
  # treated unit
  if (placebo_graphs == T) {
    
    return(list(w = w_final, 
                Y_true = Y_true_final, 
                Y_elast = Y_elast_final, 
                fit = fit_final, 
                placebo_frame = new_frame,
                placebo = placebo,
                path.plot = path.plot,
                cv_results = cv_results,
                p_stat = p_stat,
                ratio.plot = ratio.plot))
  }
  
  if (placebo_graphs == F) {
    
    return(list(w = w_final, 
                Y_true = Y_true_final, 
                Y_elast = Y_elast_final, 
                fit = fit_final, 
                placebo_frame = new_frame,
                path.plot = path.plot,
                cv_results = cv_results,
                p_stat = p_stat,
                ratio.plot = ratio.plot))
    
  }
  
  
  
}
