#' ElasticSynth
#'
#' A synthetic controls method for comparitive case studies
#' that utilizes elastic net regression -- glmnet -- to assign unit weights. 
#' Utilizes a cross-validation mechanism for selecting optimal L1 and L2 regularization parameters, as well as lambda.
#' 
#' @param PredictorMatrix a p x n matrix where each row is a seperate control variable to be included in the analysis for each unit
#' @param OutcomeMatrix a t x n matrix where each row is the outcome variable of interest at time t for each unit
#' @param treated the column index for the treated unit in both PredictorMatrix and OutcomeMatrix. For multiple treated units, see 'ElasticSynthRun'
#' @param pre a vector of row indices indicating the pre period (Ex. 1:40). This is the period over which the algorithm will optimize weights.
#' @param post a vector of row indices indicating the post period (Ex. 41:50)
#' @param lambda_grid glmnet param - a vector of values for lambda to fitted at each alpha value during cross validation. For more information, see glmnet
#' @param start_month the start month of the case study (first pre-period). This is used for plotting
#' @param end_month the end month of the case study (last post-period). This is used for plotting
#' @param time_unit The gap between observations in the case study. No default. A character string, containing one of "day", "week", "month", "quarter" or "year". This can optionally be preceded by an integer and a space, or followed by "s". Ex '3 months' or 'quarter'
#' @param penaltyMat glmnet param - (optional) a n-1 x n matrix of additional penalties placed upon glmnet for selecting certain units. Useful in the case of multiple treated units. Each column of the matrix represents the penalties to be placed on the control units. The treated unit is excluded from the column. See glmnet param penalty.factor for more details. Default is no additional penalties on any unit.
#' @param max_number_units glmnet param - The maximum number of units glmnet will select. Useful for forcing parsimonious models. Default allows all units to receive weight.
#' @param lower_limit_weights The lower limit value of weight that can be placed on any unit. Default is zero (non-negative constraint). Change if desired, but be wary of overfitting.
#' @return list containing output weights for treated unit (w_final), the actual outcome values (Y_true), the fitted outcome values (Y_elast), optimal value of lambda (lambda_opt), optimal value of alpha, (alpha_opt), a dataframe of the results of placebo test (placebo_frame), a plot of the path of the treated vs. actual unit (path.plot), and a plot of the placebo test results for the treated unit. 
#' @export

ElasticSynth <- function(PredictorMatrix, 
                         OutcomeMatrix, 
                         treated, 
                         pre, 
                         post,   
                         lambda_grid = c(seq(from = 1e-03, to = 1e-01, by = 1e-02),seq(from = 2e-01, to = 100, by = 1e-01), seq(from = 200, to = 50000, by = 100)),
                         a_grid =seq(from = 0.1, to = 0.9, by = 0.01),
                         start_month,
                         end_month,
                         time_unit,
                         penaltyMat =  matrix(rep(1, ncol(OutcomeMatrix)*(ncol(OutcomeMatrix) - 1)), 
                                              ncol = ncol(OutcomeMatrix), 
                                              nrow = ncol(OutcomeMatrix) - 1),
                         max_number_units = ncol(OutcomeMatrix),
                         lower_limit_weights = 0) 
{
  
  
  suppressMessages(library(glmnet))
  suppressMessages(library(ggplot2))
  suppressMessages(library(lubridate))
  
  #### functions
  penaltyCleanup <- function(penalty, treated, i) {
    if (i <= treated & i > 1) {
      penaltyTreated = penalty[treated - 1]
      penaltyMinusTreated = penalty[-(treated - 1)]
      penalty = c(penaltyTreated, penaltyMinusTreated)
    }
    else if (i > treated & i > 1) {
      penaltyTreated = penalty[treated]
      penaltyMinusTreated = penalty[-(treated)]
      penalty = c(penaltyTreated, penaltyMinusTreated)
    }
    return(penalty)
  }
  
  #### necessary data structure
  # Predictors for units minus treated
  X0 <- as.matrix(PredictorMatrix[,-treated])
  # Predictors for treated
  X1 <- as.matrix(PredictorMatrix[,treated, drop = F])
  # Combined
  X  <- as.matrix(cbind(X1, X0))
  # Time Periods for units minus treated
  Y0 <- as.matrix(OutcomeMatrix[,-treated])
  # Time Periods for treated
  Y1 <- as.matrix(OutcomeMatrix[,treated, drop = F])
  # Combined
  Y  <- as.matrix(cbind(Y1, Y0))
  # Pre period units minus treated
  Z0 <- as.matrix(Y0[pre,])
  # Pre period for treated
  Z1 <- as.matrix(Y1[pre,1, drop = F])
  # Combined
  Z  <- as.matrix(cbind(Z1, Z0))
  
  ##### Params
  
  # Number of Predictors
  K  <- dim(X)[1]
  # Number of units
  N  <- dim(Y)[2]
  # Total time periods
  Time <- dim(Y)[1]
  # Pre period
  T0   <- dim(Z)[1]
  # Post Period
  T1   <- Time - T0
  
  div  <- as.matrix(apply(X, 1, sd))
  X    <- X / div[,rep(1, N)]
  
  ## Month frame for plotting
  start_month <- as.Date(start_month)
  end_month <- as.Date(end_month)
  month_seq  <- seq(start_month, end_month, by = time_unit)
  month_join <- data.frame(time = 1:max(post), month = month_seq)
  
  ## Find the optimal elast
  # Iterate over i
  nlambda     <- length(lambda_grid)
  na          <- length(a_grid)
  err_alpha   <- matrix(0, nrow = na, ncol = 1)
  lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1)
  
  cat('*** Searching for Optimal Lambda & Alpha ***\n')
  penaltyMat  <- cbind(penaltyMat[,treated], penaltyMat[,-treated])
  
  for (j in 1:na) {
    a         <- a_grid[j]
    cat('a =', toString(a), '\n')
    err       <- matrix(0, nrow = N, ncol = nlambda)
    for (i in 1:N) {
      Y1      <- as.matrix(Y[,i, drop = F])
      Y0      <- as.matrix(Y[,-c(i)])
      Z1      <- as.matrix(Z[,i, drop = F])
      Z0      <- as.matrix(Z[,-c(i)])
      X1      <- as.matrix(X[,i, drop = F])
      X0      <- as.matrix(X[,-c(i)])
      Z1_tr   <- Z1
      Z0_tr   <- Z0
      Z1_te   <- as.matrix(Y1[-(1:(T0-5)),])
      Z0_te   <- as.matrix(Y0[-(1:(T0-5)),])

      # CV - Find optimal lambda and alpha across all units
      V1      <- rbind(scale(Z1_tr, scale = FALSE), scale(X1, scale = FALSE))
      V0      <- rbind(scale(Z0_tr, scale = FALSE), scale(X0, scale = FALSE))
      penalty <- penaltyCleanup(penaltyMat[,i, drop = F], treated, i)
      fit     <- glmnet(x = V0, y = V1,
                        alpha = a,
                        lambda = lambda_grid,
                        standardize = FALSE,
                        intercept = FALSE,
                        lower.limits = lower_limit_weights,
                        pmax = max_number_units,
                        penalty.factor = penalty)
      w       <- as.matrix(coef(fit, s = lambda_grid))[-1,]
      int     <- t(as.matrix(apply(Z1_tr[,rep(1, nlambda)] - Z0_tr %*% w, 2, mean)))
      e       <- Z1_te[,rep(1, nlambda)] - int[rep(1, T1 + 5),] - Z0_te %*% w
      err[i,] <- colMeans(e^2, na.rm = TRUE)
    }
    # Optimal lambda
    err       <- apply(err, 2, mean)
    ind_opt   <- which.min(err)
    err_alpha[j] <- err[ind_opt]
    lambda_opt_alpha[j] <- lambda_grid[ind_opt]
  }
  # Optimal a
  ind_opt     <- which.min(err_alpha)
  a_opt       <- a_grid[ind_opt]
  lambda_opt  <- lambda_opt_alpha[ind_opt]
  
  

  int_elast   <- matrix(0, nrow = 1, ncol = 1)
  Y_elast     <- matrix(0, nrow = Time, ncol = 1)
  Y_true      <- matrix(0, nrow = Time, ncol = 1)
  old_frame   <- data.frame()
    for (i in 1:N) {
      if (i != 1) {
        old_frame <- new_frame
      }
      Y1 <- as.matrix(Y[,i,drop = F])
      Y0 <- as.matrix(Y[,-i])
      Z1 <- as.matrix(Z[,i, drop = F])
      Z0 <- as.matrix(Z[,-i])
      X1 <- as.matrix(X[,i, drop = F])
      X0 <- as.matrix(X[,-i])
      V1 <- rbind(scale(Z1, scale = FALSE), scale(X1, scale = FALSE))
      V0 <- rbind(scale(Z0, scale = FALSE), scale(X0, scale = FALSE))
      penalty = penaltyCleanup(penaltyMat[,i, drop = F], treated, i)
      
      # Fit elast
      fit_final   <- glmnet(x = V0, y = V1,
                            alpha = a_opt,
                            lambda = lambda_grid,
                            standardize = FALSE,
                            intercept = FALSE,
                            pmax = max_number_units,
                            lower.limits = lower_limit_weights,
                            penalty.factor = penalty)
      w           <- as.matrix(coef(fit_final, s = lambda_opt))[-1,]
      if (i == 1) {
        w_final     <- w
        int_elast   <- as.matrix(apply(Z1 - Z0 %*% w_final, 2, mean))
        Y_elast     <- int_elast[rep(1, Time),] + Y0 %*% w_final
        Y_elast_final <- Y_elast
        Y_true      <- Y1[c(pre,post)]
        Y_true_final <- Y_true
        gaps        <- Y_true[c(pre,post)] - Y_elast[c(pre,post)]
        rmse        <- sqrt(sum(abs(gaps[pre]))/length(pre))
        
        
        plotFrame   <- data.frame(time = c(pre,post), Y_true = Y_true, Y_elast = Y_elast, gaps = gaps)
        plotFrame   <- merge(plotFrame, month_join, by = 'time', all.x = TRUE)
        min_post_month <- plotFrame$month[plotFrame$time == min(post)]
        max_post_month <- plotFrame$month[plotFrame$time == max(post)]
        # 
        path.plot   <- ggplot(plotFrame) +
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
      
      int_elast   <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
      Y_elast     <- int_elast[rep(1, Time),] + Y0 %*% w
      Y_true      <- Y1
      gaps        <- Y_true[c(pre,post)] - Y_elast[c(pre,post)]
      rmse        <- sqrt(sum(abs(gaps[pre]))/length(pre))
      new_frame   <- data.frame(gaps = unlist(gaps), time = c(pre,post), unit_type = ifelse(i == 1, 'Treated Unit', 'Control Unit Distribution'), unit = colnames(Y)[i])
      new_frame   <- rbind(old_frame, new_frame)
      
    }
    new_frame     <- merge(new_frame, month_join, by = 'time', all.x = TRUE)
    min_post_month <- unique(new_frame$month[new_frame$time == min(post)])
    max_post_month <- unique(new_frame$month[new_frame$time == max(post)])
    far_x_axis_month <- max_post_month %m+% months(2)
    near_x_axis_month <- max_post_month %m-% months(3*10)
    
    
    placebo <- ggplot() +
      geom_line(data = subset(new_frame, unit_type == 'Control Unit Distribution'), aes(x = month, y = gaps, group = unit),colour = 'grey', size = 1) +
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
                lambda_opt = lambda_opt, 
                alpha_opt = a_opt,
                placebo_frame = new_frame,
                placebo = placebo,
                path.plot = path.plot))
    
    
}



