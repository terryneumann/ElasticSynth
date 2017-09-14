
#' ElasticSynthRun
#'
#' An implementation of ElasticSynth that can run on multiple treated units, through multiple pre and post periods. No output object is returned, 
#' but graphical results and output weights are stored in output directory.
#' 
#' @param PredictorMatrix a p x n matrix where each row is a seperate control variable to be included in the analysis for each unit
#' @param OutcomeMatrix a t x n matrix where each row is the outcome variable of interest at time t for each unit
#' @param TreatedUnitList a vector/list of column indices of the treated unit 
#' @param PrePeriodList a list of lists of the pre periods for the treated units. Ex. list(list(1:40), list(1:39), list(2:39))
#' @param PostPeriodList a list of lists of the post periods for the treated units. Ex list(list(41:50), list(40:50), list(40:50))
#' @param max_number_units glmnet param - The maximum number of units glmnet will select. Useful for forcing parsimonious models. Default allows all units to receive weight.
#' @param lower_limit_weights glmnet param - The lower limit value of weight that can be placed on any unit. Default is zero (non-negative constraint). Change if desired, but be wary of overfitting.
#' @param PenaltyTermsMat glmnet param - (optional) a n-1 x n matrix of additional penalties placed upon glmnet for selecting certain units. Useful in the case of multiple treated units. Each column of the matrix represents the penalties to be placed on the control units. The treated unit is excluded from the column. See glmnet param penalty.factor for more details. Default is no additional penalties on any unit.
#' @param start_month the start month of the case study (first pre-period). This is used for plotting
#' @param end_month the end month of the case study (last post-period). This is used for plotting
#' @param time_unit The gap between observations in the case study. No default. A character string, containing one of "day", "week", "month", "quarter" or "year". This can optionally be preceded by an integer and a space, or followed by "s". Ex '3 months' or 'quarter'
#' @param storage_folder The directory where graph and weight outputs should be stored
#' @return Doesn't return an object, but will save placebo test results and true vs. fitted graphs for each treated unit into storage_folder param

ElasticSynthRun <- function(PredictorMatrix,
                            OutcomeMatrix,
                            TreatedUnitList,
                            PrePeriodList,
                            PostPeriodList,
                            max_number_units,
                            PenaltyTermsMat = matrix(rep(1, ncol(OutcomeMatrix)*(ncol(OutcomeMatrix) - 1)), 
                                                     ncol = ncol(OutcomeMatrix), 
                                                     nrow = ncol(OutcomeMatrix)),
                            time_unit,
                            start_month,
                            end_month,
                            lower_limit_weights,
                            storage_folder) 

{
  
  suppressMessages(library(ggplot2))
  
  
  n          <- length(TreatedUnitList)
  weights    <- rep(list(list()), n)
  
  for (i in 1:n) {
    cat(paste('
              ##########\n
              ##########\n
              ##########\n
              Treated Unit', TreatedUnitList[i]))
    elasticSynth.out <- ElasticSynth(PredictorMatrix = PredictorMatrix,
                                     OutcomeMatrix = OutcomeMatrix,
                                     treated = TreatedUnitList[i],
                                     pre = unlist(PrePeriodList[[i]]),
                                     post = unlist(PostPeriodList[[i]]),
                                     max_number_units = max_number_units,
                                     time_unit = time_unit,
                                     lambda_grid = c(seq(from = 1e-03, to = 1e-01, by = 1e-02),seq(from = 2e-01, to = 100, by = 1e-01), seq(from = 200, to = 50000, by = 100)),
                                     a_grid =seq(from = 0.1, to = 0.9, by = 0.01),
                                     start_month = start_month,
                                     end_month = end_month,
                                     lower_limit_weights = lower_limit_weights,
                                     penaltyMat = PenaltyTermsMat
    )
    
    
    placebo <- elasticSynth.out$placebo
    path.plot <- elasticSynth.out$path.plot 
    
    
    ggsave(path.plot, filename = paste(storage_folder, '/path', TreatedUnitList[i], '.png', sep = ''), width = 10, height = 5)
    ggsave(placebo, filename = paste(storage_folder, '/placebo', TreatedUnitList[i], '.png', sep = ''), width = 10, height = 5)
    
    
    weights[[i]] <- elasticSynth.out$w
    
    
  }
  
  
  save(weights, file = paste(storage_folder, '/SynthWeights.Rdata', sep = ''))
  
}
