
#' ElasticSynth - a class for running synthetic controls with elastic net regularization
#'
#'
#' Model is initiated with a list of parameters most which include treated_units and donor_units, two dataframes that 
#' contain time periods, units, and counts of outcomes
#' 
#' @param treated_units Treated units data.frame
#' @param donor_units Donor units data.frame
#' @param measure_vars Outcome measures upon which to jointly fit a synthetic control
#' @param measure_col Column name for measure_vars
#' @param unit_col Column name for unit information
#' @param time_col Column name for time information (in form 1:max(time))
#' @param pre_list A list of lists or vectors representing pre periods across all outcomes
#' @param post_list A list of lists or vectors representing post periods across all outcomes
#' @param start_cv_list A list indicating where to start CV procedure in time
#' @param end_cv_list A list indicating where to end CV procedure in time
#' @param time_unit_list A list indicating the time between time periods
#' @param alphas vector of alphas to use for cross validation - for more information see glment
#' @importFrom R6 R6Class
#' @export
ElasticSynth <- R6::R6Class(
  class = FALSE,
  private = list(
    cv_results = NULL,
    weights = NULL,
    series_frame = NULL,
    type = NULL
  ),
  public = list(
    treated_units = NULL,
    donor_units = NULL,
    measure_vars = NULL,
    measure_col = NULL,
    unit_col = NULL,
    time_col = NULL,
    value_col = NULL,
    pre_list = NULL,
    post_list = NULL,
    start_cv_list = NULL,
    end_cv_list = NULL,
    cv_step_list = NULL,
    start_month_list = NULL,
    end_month_list = NULL,
    time_unit_list = NULL,
    alphas = NULL,
    lambdas = NULL,
    lower_limit_weights = NULL,
    upper_limit_weights = NULL,
    placebo = NULL,
    verbose = NULL,
    treated = NULL,
    max_pre_month = NULL,
    end_month = NULL,
    start_month = NULL,
    class_name = 'ElasticSynth',
    initialize = function(treated_units,
                         donor_units,
                         measure_vars,
                         measure_col,
                         unit_col,
                         time_col,
                         value_col,
                         pre_list,
                         post_list,
                         start_cv_list,
                         end_cv_list,
                         cv_step_list,
                         start_month_list,
                         end_month_list,
                         time_unit_list,
                         alphas = NULL,
                         lambdas = NULL,
                         lower_limit_weights,
                         upper_limit_weights,
                         placebo,
                         verbose,
                         max_pre_month,
                         end_month,
                         start_month
    )
    {
      ### load packages
      suppressMessages(library(glmnet))
      suppressMessages(library(ggplot2))
      suppressMessages(library(lubridate))
      suppressMessages(library(dplyr))
      suppressMessages(library(tidyr))
      suppressMessages(library(Metrics))
      suppressMessages(library(reshape2))
      suppressMessages(library(data.table))
      suppressMessages(library(ggplot2))
      
      
      if (is.null(alphas)) {
        alphas <- seq(0, 1, by = 0.05)
      }
      
      if (is.null(lambdas)) {
        lambdas <-
          c(seq(1, 0.02, by = -0.1),
            seq(0.1, 0.0001, by = -.0005),
            0.0001,
            0.00005)
      }
      
      if (is.null(lower_limit_weights)) {
        lower_limit_weights <- 0
      }
      
      if (is.null(upper_limit_weights)) {
        upper_limit_weights <- 1
      }
      
      
      ### Check types
      stopifnot(is.data.frame(treated_units))
      stopifnot(is.data.frame(donor_units))
      stopifnot(is.character(measure_vars))
      stopifnot(is.character(measure_col))
      stopifnot(is.character(unit_col))
      stopifnot(is.character(time_col))
      stopifnot(is.character(value_col))
      stopifnot(is.character(start_month_list))
      stopifnot(is.character(end_month_list))
      stopifnot(is.character(time_unit_list))
      stopifnot(is.list(pre_list))
      stopifnot(is.list(post_list))
      stopifnot(is.numeric(start_cv_list))
      stopifnot(is.numeric(end_cv_list))
      stopifnot(is.numeric(cv_step_list))
      stopifnot(is.numeric(alphas))
      stopifnot(is.numeric(lambdas))
      stopifnot(is.numeric(lower_limit_weights))
      stopifnot(is.numeric(upper_limit_weights))
      stopifnot(lower_limit_weights < upper_limit_weights)
      # -- clean this up
      stopifnot(lubridate::is.Date(max_pre_month))
      stopifnot(lubridate::is.Date(start_month))
      stopifnot(lubridate::is.Date(end_month))
      # --
      stopifnot(is.logical(placebo))
      stopifnot(is.logical(verbose))
      
      ### Assign
      self$treated <- 1
      self$alphas <- alphas
      self$lambdas <- lambdas
      self$treated_units <- treated_units
      self$donor_units <- donor_units
      self$measure_vars <- measure_vars
      self$measure_col <- measure_col
      self$unit_col <- unit_col
      self$time_col <- time_col
      self$value_col <- value_col
      self$pre_list <- pre_list
      self$post_list <- post_list
      self$start_month_list <- start_month_list
      self$end_month_list <- end_month_list
      self$start_cv_list <- start_cv_list
      self$end_cv_list <- end_cv_list
      self$time_unit_list <- time_unit_list
      self$cv_step_list <- cv_step_list
      self$max_pre_month <- max_pre_month
      self$start_month <- start_month
      self$end_month <- end_month
      self$lower_limit_weights <- lower_limit_weights
      self$upper_limit_weights <- upper_limit_weights
      self$placebo <- placebo
      self$verbose <- verbose
    },
    
    cv_treated = function() {
      
      treated_units = setDT(self$treated_units)
      donor_units = setDT(self$donor_units)
      raw_results = data.frame()
      for (m in 1:length(self$measure_vars)) {
        cat(
          paste(
            '\n\n\n\n\n\n\nCross Validating for measure:',
            self$measure_vars[m],
            '\n\n\n\n\n'
          )
        )
        pre = self$pre_list[[m]]
        post = self$post_list[[m]]
        
        treated_wide_measure = self$long_to_wide(self$treated_units[base::get(self$measure_col) == self$measure_vars[m], ],
                                            self$time_col,
                                            self$unit_col,
                                            self$value_col)
        donor_wide_measure = self$long_to_wide(self$donor_units[base::get(self$measure_col) == self$measure_vars[m], ],
                                          self$time_col,
                                          self$unit_col,
                                          self$value_col)
        # generate a list of the hold out periods based on Hyndman CV structure
        folds_ix = self$generate_hyndman_folds(
          start_cv = self$start_cv_list[m],
          end_cv = self$end_cv_list[m],
          cv_step = self$cv_step_list[m]
        )
        ####  data structure
        OutcomeMatrix = as.matrix(cbind(treated_wide_measure, donor_wide_measure))
        # Time Periods for units minus treated
        Y0 = as.matrix(OutcomeMatrix[, -self$treated, drop = F])
        # Time Periods for treated
        Y1 = as.matrix(OutcomeMatrix[, self$treated, drop = F])
        # Combined
        Y  = as.matrix(cbind(Y1, Y0))
        # Pre period units minus treated
        Z0 = as.matrix(Y0[pre, ])
        # Pre period for treated
        Z1 = as.matrix(Y1[pre, 1, drop = F])
        # Combined
        Z  = as.matrix(cbind(Z1, Z0))
        
        
        # Number of units
        if (self$placebo == T) {
          N  = dim(Y)[2]
        }
        
        if (self$placebo == F) {
          N = ncol(treated_wide_measure)
        }
        
        ## Month frame for plotting
        na  = length(self$alphas)
        
        cat('\n\n\n\n******* Unit Cross Validation *******\n\n\n\n')
        
        units_for_same_measure = data.frame()
        for (i in 1:N) {
          if (self$verbose == T) {
            cat('*** Unit', colnames(Z[, i, drop = F]), '***\n\n')
          }
          
          Z1       = as.matrix(Z[, i, drop = F])
          Z0       = as.matrix(Z[, -c(1:ncol(treated_wide_measure), i), drop = F])
          cv_grid = data.frame()
          for (r in 1:length(folds_ix)) {
            # hyndman cv structure
            Z0_train = Z0[1:(min(folds_ix[[r]]) - 1), ]
            Z1_train = Z1[1:(min(folds_ix[[r]]) - 1), ]
            Z0_test  = Z0[min(folds_ix[[r]]):max(folds_ix[[r]]), ]
            Z1_test  = Z1[min(folds_ix[[r]]):max(folds_ix[[r]]), ]
            alpha_grid = data.frame()
            for (j in 1:na) {
              a        = self$alphas[j]
              fit      = glmnet(
                x = Z0_train,
                y = Z1_train,
                family = 'gaussian',
                alpha = a,
                lambda = self$lambdas,
                standardize = FALSE,
                intercept = FALSE,
                lower.limits = self$lower_limit_weights,
                upper.limits = self$upper_limit_weights
              )
              
              w = as.matrix(coef(fit, s = self$lambdas))[-1, ]
              int_elast = as.matrix(apply(Z1_train - Z0_train %*% w, 2, mean))
              y_pred  = int_elast[rep(1, nrow(Z0_test)), ]  + as.matrix(Z0_test) %*% w
              y_true  = Z1_test
              # Chose MAPE because it is unitless across multiple outcomes
              mape_lambda = as.numeric(apply(y_pred, 2, function(y_pred)
                Metrics::mape(y_true, y_pred)))
              alpha_grid_tmp = data.frame(
                fold = paste0('Fold', rep(r, length(self$lambdas))),
                alpha = rep(a, length(self$lambdas)),
                lambda = self$lambdas,
                mape = mape_lambda
              )
              alpha_grid = rbind(alpha_grid, alpha_grid_tmp)
            }
            alpha_grid$measure = self$measure_vars[m]
            alpha_grid$unit    = colnames(Z)[i]
            cv_grid = rbind(cv_grid, alpha_grid)
          }
          names(cv_grid)[which(names(cv_grid) == self$measure_col)] = 'measure'
          wide_cv_grid = reshape2::dcast(cv_grid, lambda + alpha + measure + unit ~ fold, value.var = 'mape')
          names(wide_cv_grid)[1:4] = c('lambda', 'alpha', 'measure', 'unit')
          # results for all units across same measure
          units_for_same_measure = rbind(units_for_same_measure, wide_cv_grid)
        }
        # results for all units across all measures
        raw_results = bind_rows(units_for_same_measure, raw_results)
      }
      
      summary_results = raw_results %>%
        mutate(row_means = rowMeans(dplyr::select(., starts_with('Fold')), na.rm = T)) %>%
        group_by(lambda, alpha, unit) %>%
        summarise(cv_score = mean(row_means, na.rm = T)) %>%
        ungroup() %>%
        group_by(unit) %>%
        summarise(alpha = alpha[which.min(cv_score)],
                  lambda = lambda[which.min(cv_score)],
                  cv_score = min(cv_score)) %>%
        ungroup()
      
      summary_results$lambda_ix = 0
      for (u in 1:nrow(summary_results)) {
        summary_results$lambda_ix[u] = which(self$lambdas == summary_results$lambda[u])
      }
      summary_results = summary_results %>%
        arrange(factor(unit, levels = colnames(Y)))
      return(summary_results)
    },
    cv_untreated = function() {
      setDT(self$donor_units)
      setDT(self$treated_units)
      results = data.frame()
      for (m in 1:length(self$measure_vars)) {
        measure_results <- data.frame()
        cat(
          paste(
            '\n\n\n\n\n\n\nCross Validating for measure:',
            self$measure_vars[m],
            '\n\n\n\n\n'
          )
        )
        pre = self$pre_list[[m]]
        post = self$post_list[[m]]
        
        donor_wide_measure = self$long_to_wide(self$donor_units[base::get(self$measure_col) == self$measure_vars[m], ],
                                               self$time_col,outcome_join_finalget
                                               self$unit_col,
                                               self$value_col)
        # generate a list of the hold out periods based on Hyndman CV structure
        folds_ix = self$generate_hyndman_folds(
          start_cv = self$start_cv_list[m],
          end_cv = self$end_cv_list[m],
          cv_step = self$cv_step_list[m]
        )
        ####  data structure
        Y = as.matrix(donor_wide_measure)
        
        # Number of units
        N  = dim(Y)[2]
        
        
        
        # Total time periods
        Time = length(c(pre,post))
        # Pre period
        T0   = length(pre)
        # Post Period
        T1   = length(post)
        na  = length(self$alphas)
        
        
        cat('\n\n\n\n******* Cross Validation *******\n\n\n\n')
        fold_grid = data.frame()
        for (r in 1:length(folds_ix)) {
          cat(paste0('***************** fold ', r, '\n'))
          cv_grid = data.frame()
          for (i in 1:N) {
            # unitName = sub('District ', '', colnames(Y)[i])
            
            
            Z1       = as.matrix(Y[,i, drop = F])
            Z0       = as.matrix(Y[,-i,drop = F])   
            Z0_train = Z0[1:(min(folds_ix[[r]]) - 1),,drop=F] 
            Z1_train = Z1[1:(min(folds_ix[[r]]) - 1),,drop=F]
            Z0_test  = Z0[min(folds_ix[[r]]):max(folds_ix[[r]]),,drop=F]
            Z1_test  = Z1[min(folds_ix[[r]]):max(folds_ix[[r]]),,drop=F]
            
            if (self$verbose == T) {
              cat('*** Unit', colnames(Y[,i,drop = F]), '***\n\n')
            }
            
            alpha_grid = data.frame()
            # DI cv structure
            for (j in 1:na) {
              a        = self$alphas[j]
              fit      = glmnet(x = Z0_train, 
                                y = Z1_train,
                                family = 'gaussian',
                                alpha = a,
                                lambda = self$lambdas,
                                standardize = FALSE,
                                intercept = FALSE,
                                lower.limits = self$lower_limit_weights,
                                upper.limits = self$upper_limit_weights)
              y_pred  = predict(fit, newx = Z0_test, type = 'response', s = fit$lambda)
              y_true  = Z1_test
              mape_lambda = as.numeric(apply(y_pred, 2, function(y_pred)
                Metrics::mape(y_true, y_pred)))
              lambda_ix  = which(mape_lambda == min(mape_lambda))
              alpha_grid_tmp = data.frame(unit = colnames(Z1_train), 
                                          fold = rep(r, length(self$lambdas)), 
                                          alpha = rep(a, length(self$lambdas)), 
                                          lambda = self$lambdas, 
                                          mape = mape_lambda,
                                          measure = self$measure_vars[m])
              alpha_grid = rbind(alpha_grid, alpha_grid_tmp)
            }
            cv_grid = rbind(cv_grid, alpha_grid)
            rm(alpha_grid)
          }
          fold_grid = rbind(cv_grid, fold_grid)
          rm(cv_grid)
        }
        fold_grid[fold_grid == Inf] = NA
        fold_grid = fold_grid %>%
          mutate(fold= paste0('fold_', fold)) %>%
          group_by(alpha, lambda, measure, fold) %>%
          summarise(mape = mean(mape, na.rm = T))
        
        final_grid = reshape2::dcast(fold_grid, lambda + alpha ~ fold + measure, value.var = 'mape')
        rm(fold_grid)
        # results for all units across all measures
        if (m==1) {
          results = final_grid
        } else {
          results = results %>%
            left_join(final_grid)
          
        }
      }
      
      results = results %>%
        mutate(rowmeans = rowMeans(results[,3:ncol(results)]))
      min_row = which(results$rowmeans == min(results$rowmeans))
      choice_alpha = results[min_row,'alpha']
      choice_lambda = results[min_row, 'lambda']
      choice_lambda_ix = which(self$lambdas == choice_lambda)[1]
      cat(paste('choice alpha', choice_alpha, 'choice lambda', choice_lambda, 'choice lambda ix', choice_lambda_ix))
      units = c(unique(as.character(self$treated_units[[self$unit_col]])), 
                unique(as.character(self$donor_units[[self$unit_col]])))
      cat(paste('units', units))
      summary_results = expand.grid(unit = units, 
                                    alpha = choice_alpha, 
                                    lambda = choice_lambda,
                                    lambda_ix = choice_lambda_ix,
                                    mape = results[min_row, 'rowmeans'])
      return(summary_results)
    },
    cv = function(type) {
      private$type <- type
      if (private$type == 'individual') {
        results <- self$cv_treated()
      } else if (private$type == 'grouped') {
        results <- self$cv_untreated()
      }
      return(results)
    },
    generate_weights = function(cv_results) {
      
      ## Pre periods for stacked outcomes
      treated_units_wide = self$long_to_wide_weights(dt = self$treated_units, 
                                                self$pre_list, 
                                                self$time_col, 
                                                self$unit_col, 
                                                self$measure_col, 
                                                self$value_col, 
                                                self$measure_vars)
      
      donor_units_wide = self$long_to_wide_weights(dt = self$donor_units, 
                                              self$pre_list, 
                                              self$time_col, 
                                              self$unit_col, 
                                              self$measure_col, 
                                              self$value_col, 
                                              self$measure_vars)
      
      ####  data structure
      OutcomeMatrix = as.matrix(cbind(treated_units_wide, donor_units_wide))
      # Time Periods for units minus treated
      Y0 = as.matrix(OutcomeMatrix[,-self$treated])
      # Time Periods for treated
      Y1 = as.matrix(OutcomeMatrix[,self$treated, drop = F])
      # Combined
      Y  = as.matrix(cbind(Y1, Y0))
      ##### Params
      n_tr = dim(self$treated_units)[2]
      # Number of units
      private$cv_results <- cv_results
      
      private$cv_results = private$cv_results %>%
        arrange(factor(unit, levels = colnames(Y)))
      
      
      if (self$placebo) {
        N = dim(Y)[2]
      }
      else {
        N = ncol(treated_units_wide)
      }
      
      w_final = list(list())

      for (i in 1:N) {
        
        
        Z1 = as.matrix(Y[,i, drop = F])
        Z0 = as.matrix(Y[,-c(1:ncol(treated_units_wide), i),drop = F])
        Z = cbind(Z1, Z0)
        row_means = apply(Z, 1, mean)
        row_max   = apply(Z, 1, max)
        
        # Still experimental --- if there are very different levels across outcomes, glmnet will over fit to the outcome with the highest levels
        # ---------------------- as it is trying to minimize overall MSE. Therefore, if there are multiple outcomes, need some way to normalize data. 
        # ---------------------- Here, I assigne observation weights of 1/row mean to normalize. Still up for debate what the best normailization approach is.
        
        if (length(self$measure_vars) == 1) {
          obs.weights = rep(1, nrow(Z0))
        }

        else {
          obs.weights = 1/row_means
        }

        
        # 
        # Fit elast
        fit_final   = glmnet(x = Z0, 
                             y = Z1,
                             family = 'gaussian',
                             alpha = private$cv_results$alpha[i],
                             lambda = self$lambdas,
                             standardize = FALSE,
                             intercept = FALSE,
                             lower.limits = self$lower_limit_weights,
                             upper.limits = self$upper_limit_weights,
                             weights = obs.weights)
        
        w_final[[i]] = as.matrix(coef(fit_final, s = self$lambdas))[-1,private$cv_results$lambda_ix[i]]
        
      }
      
      if (self$placebo == T) {
        names(w_final) = private$cv_results$unit
      }
      else {
        names(w_final) = private$cv_results$unit[1:ncol(treated_units_wide)]
      }
      w_final
    },
    fit = function(weights, cv_results) {

      setDT(self$treated_units)
      setDT(self$donor_units)
      private$cv_results = cv_results
      private$cv_results$fitted_error_pre = 0
      private$cv_results$fitted_mape_pre  = 0
      private$cv_results$error_post = 0
      private$cv_results$mape_post = 0
      private$weights = weights
      result_cv_frame = data.frame()
      result_series_frame = data.frame()
      
      for (j in 1:length(self$measure_vars)) {
        
        
        start_month = as.Date(self$start_month_list[j])
        end_month = as.Date(self$end_month_list[j])
        month_seq  = seq(start_month, end_month, by = self$time_unit_list[j])
        month_join = data.frame(time = min(self$pre_list[[j]]):max(self$post_list[[j]]), month = month_seq)
        treated_units_wide = self$long_to_wide(self$treated_units[base::get(self$measure_col) == self$measure_vars[j],], 
                                          self$time_col, 
                                          self$unit_col, 
                                          self$value_col)
        donor_units_wide = self$long_to_wide(self$donor_units[base::get(self$measure_col) == self$measure_vars[j],], 
                                        self$time_col, 
                                        self$unit_col, 
                                        self$value_col)
        
        OutcomeMatrix = as.matrix(cbind(treated_units_wide, donor_units_wide))
        Y0 = as.matrix(OutcomeMatrix[,-self$treated])
        Y1 = as.matrix(OutcomeMatrix[,self$treated, drop = F])
        Y  = as.matrix(cbind(Y1, Y0))
        
        private$cv_results = private$cv_results %>%
          arrange(factor(unit, levels = colnames(Y)))
        
        ##### Params
        n_tr = dim(self$treated_units)[2]
        # Number of units
        
        # Pre period
        T0   = length(self$pre_list[[j]])
        # Post Period
        T1   = length(self$post_list[[j]])
        
        if (self$placebo) {
          N = dim(OutcomeMatrix)[2]
        }
        else {
          N = ncol(treated_units_wide)
        }
        
        series_frame = data.frame()
        
        for (i in 1:N) {
          Y0_res = as.matrix(Y[,-c(1:ncol(treated_units_wide), i)])
          # Time Periods for treated
          Y1_res = as.matrix(Y[,i, drop = F])
          Y = cbind(Y1, Y0)
          Z0 = Y0_res[pre_list[[j]],,drop=F]
          Z1 = Y1_res[pre_list[[j]],,drop=F]
          
          int_elast        = as.matrix(apply(Z1 - as.matrix(Z0[self$pre_list[[j]],]) %*% as.matrix(private$weights[[i]]), 2, mean))
          Y_elast          = int_elast[rep(1, max(self$post_list[[j]])),] + as.matrix(Y0_res) %*% private$weights[[i]]
          Y_true           = Y1_res
          gaps             = Y_true - Y_elast
          private$cv_results$fitted_error_pre[i] = Metrics::mse(Y_true[self$pre_list[[j]]], Y_elast[self$pre_list[[j]]])
          private$cv_results$fitted_mape_pre[i]  = Metrics::mape(Y_true[self$pre_list[[j]]], Y_elast[self$pre_list[[j]]])
          private$cv_results$error_post[i] = Metrics::mse(Y_true[self$post_list[[j]]], Y_elast[self$post_list[[j]]])
          private$cv_results$mape_post[i] = Metrics::mape(Y_true[self$post_list[[j]]], Y_elast[self$post_list[[j]]])
          

          series_frame   = rbind(series_frame, 
                                 data.frame(gaps = unname(unlist(gaps)), 
                                            y_true = Y_true[c(min(self$pre_list[[j]]):max(self$pre_list[[j]]),self$post_list[[j]])],
                                            y_elast = Y_elast[c(min(self$pre_list[[j]]):max(self$pre_list[[j]]),self$post_list[[j]])],
                                            time = c(min(self$pre_list[[j]]):max(self$pre_list[[j]]),self$post_list[[j]]), 
                                            unit_type = ifelse(i %in% 1:ncol(treated_units_wide), 'Treated Unit', 'Control Unit Distribution'), 
                                            unit = colnames(Y)[i]))
          
        }
        
        series_frame     = merge(series_frame, month_join, by = 'time', all.x = TRUE)
        private$cv_results$placebo_test_p_value = 1 - rank(private$cv_results$error_post)/nrow(private$cv_results)
        private$cv_results$adh_test_p_value = 1 - rank(private$cv_results$mape_post/private$cv_results$fitted_mape_pre)/nrow(private$cv_results)
        private$cv_results$series = self$measure_vars[j]
        series_frame$series = self$measure_vars[j]
        
        result_cv_frame = rbind(result_cv_frame, private$cv_results)
        result_series_frame = rbind(result_series_frame, series_frame)
        
      }
      
      return(list(stats = result_cv_frame,
                  series_frame = result_series_frame))
      
    },
    plot = function(series_frame) {
      private$series_frame <- series_frame
      treated_unit_names <- sort(unique(subset(private$series_frame, unit_type == 'Treated Unit')$unit))
      path_plots <- list()
      placebo_plots <- list()
      for (i in 1:length(treated_unit_names)) {
        path_plot = ggplot(subset(private$series_frame, unit == treated_unit_names[i])) +
          geom_line(aes(x = as.Date(month), y = y_true), colour = 'black', size = 2) +
          geom_line(aes(x = as.Date(month), y = y_elast), colour = 'indianred3', linetype = 'F1', size = 2, alpha = 0.7) +
          geom_vline(xintercept = as.numeric(as.Date(max_pre_month)), linetype = 'dotted', size = 1.4, alpha = 0.8) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                text = element_text(size = 10)) +
          xlab('Time Period') +
          ylab('Outcome') +  
          xlim(c(as.Date(start_month), as.Date(end_month) + 60)) +
          ggtitle(paste('Actual (Black) vs Synthetic for', treated_unit_names[i], sep = ' ')) +
          facet_wrap(~series, scales = 'free_y', nrow = 3)
        
        path_plots[[i]] <- path_plot
        
        placebo_plot = ggplot(private$series_frame) +
          geom_line(data = subset(private$series_frame, unit_type == 'Control Unit Distribution'), aes(x = month, y = gaps, group = unit),colour = 'grey', size = 0.5, alpha = 0.55) +
          geom_line(data = subset(private$series_frame, unit == treated_unit_names[i]), aes(x = month, y = gaps), size = 2, colour = 'black') +
          geom_hline(aes(yintercept = 0)) +  facet_wrap(~series, scales = 'free_y', nrow = 3) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"), 
                axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                text = element_text(size = 10)) +
          geom_vline(xintercept = as.numeric(as.Date(self$max_pre_month)), linetype = 'dotted', size = 1.4, alpha = 0.8)
        placebo_plots[[i]] <- placebo_plot
      }
      
      names(path_plots) = treated_unit_names
      names(placebo_plots) = treated_unit_names
      return(list(path_plots = path_plots, placebo_plots = placebo_plots))
    },
    generate_hyndman_folds = function(start_cv, end_cv, cv_step) {
      number_folds <-  end_cv - (start_cv - 1) - (cv_step - 1)
      folds_ix <- list()
      for (i in 1:number_folds) {
        folds_ix[[i]] <- (start_cv + (i - 1)):(start_cv + (i - 1) + (cv_step - 1))
      }
      folds_ix
    },
    long_to_wide = function(dt, time_col, unit_col, value_col) {
      dt_wide = reshape2::dcast(dt, get(time_col) ~ get(unit_col), value.var = value_col)
      dt_wide[,1] = NULL
      dt_wide[is.na(dt_wide)] = 0
      dt_wide[dt_wide == Inf] = 0
      dt_wide
    },
    long_to_wide_weights = function(dt, pre_list, time_col, unit_col, measure_col, value_col, fitted_vars) {
      library(data.table)
      setDT(dt)
      subset_frame = data.frame()
      for (i in 1:length(fitted_vars)) {
        tmp = data.frame(measure = fitted_vars[i], pre = pre_list[[i]])
        subset_frame = rbind(tmp, subset_frame)
      }
      dt = merge(dt, subset_frame, by.x = c(measure_col, time_col),by.y = c('measure', 'pre'), all.y = T) 
      dt_wide = reshape2::dcast(dt, get(time_col) + get(measure_col) ~ get(unit_col), value.var = value_col) %>% setDT()
      names(dt_wide)[1:2] = c(time_col, measure_col)
      dt_wide[,eval(measure_col) := factor(get(measure_col), levels = fitted_vars)]
      dt_wide = dt_wide %>% arrange(get(measure_col), get(time_col))
      dt_wide[,1:2] = NULL
      dt_wide[is.na(dt_wide)] = 0
      dt_wide[dt_wide == Inf] = 0
      dt_wide
    }
  )
)
    
 
 

 