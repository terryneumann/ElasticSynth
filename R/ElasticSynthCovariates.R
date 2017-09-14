#' ElasticSynthCovariates
#'
#' 
#' A function that uses the weights generated from ElasticSynth or ElasticSynthRun to explore covariates. In other words, this function
#' asks "in a synthetic unit similar to unit x in terms of outcome y, how does the synthetic unit's variable z compare to the actual unit?"
#' By finding large differences between the synthetic and actual unit, you may find what is driving the change in the intervention.
#' 
#' @param treatedUnits a list of column indices of treated units in the covariate matrix
#' @param weights either a list or a list of lists of the weights generated from ElasticSynth or ElasticSynthRun
#' @param pre a vector of pre periods. Since no optimization is happening in this step, enter the minimum pre period and the maximum pre period
#' @param post a vector of post periods. Since no optimization is happening in this step, enter the minimum post period and the maximum post period
#' @param start_month the start month of the case study (first pre-period). This is used for plotting
#' @param end_month the end month of the case study (last post-period). This is used for plotting
#' @param time_unit The gap between observations in the case study. No default. A character string, containing one of "day", "week", "month", "quarter" or "year". This can optionally be preceded by an integer and a space, or followed by "s". Ex '3 months' or 'quarter'
#' @param storage_folder The directory where graph outputs should be stored
#' @return Doesn't return an object, but will output path plot and gap plot to a directory
#' @export


ElasticSynthCovariates <- function(treatedUnits,
                                   weights,
                                   covariate,
                                   pre, 
                                   post, 
                                   start_month, 
                                   end_month,
                                   time_unit,
                                   storage_folder) 
  
{
  
  suppressMessages(library(ggplot2))
  suppressMessages(library(lubridate))
  
  
  
  n = length(treatedUnits)
  
  for (i in 1:n) {
    
    actual = covariate[,treatedUnits[i]]
    controls = covariate[,-treatedUnits[i]]
    
    synthCovariate = rowSums(as.matrix(controls) %*% diag(weights[[i]]))
    
    assign(paste('output', i, sep = ''), data.frame(month = seq(as.Date(start_month), as.Date(end_month), by = time_unit),
                                                    time = 1:max(post),
                                                    actual = actual, 
                                                    synthCovariate = synthCovariate,
                                                    gaps = actual - synthCovariate))
    
    
    min_post_month = get(paste('output', i, sep = ''))$month[get(paste('output', i, sep = ''))$time == min(post)]
    max_post_month = get(paste('output', i, sep = ''))$month[get(paste('output', i, sep = ''))$time == max(post)]
    far_x_axis_month = max_post_month %m+% months(2)
    near_x_axis_month = max_post_month %m-% months(3*10)
    
    
    
    path.plot <-     ggplot(get(paste('output', i, sep = ''))) +
      geom_line(aes(x = as.Date(month), y = actual), colour = 'black', size = 2) +
      geom_line(aes(x = as.Date(month), y = synthCovariate), colour = 'indianred3', linetype = 'longdash', size = 2) +
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
      xlim(c(as.Date(start_month), as.Date(far_x_axis_month)))
    
    
    gaps.plot <-  ggplot(get(paste('output', i, sep = ''))) +
      geom_line(aes(x = month, y = gaps), size = 2, colour = 'black') +
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
    
    
    
    
    ggsave(path.plot, filename = paste(storage_folder, '/path.plot', treatedUnits[i], '.png', sep = ''), width = 10, height = 5)     
    ggsave(gaps.plot, filename = paste(storage_folder, '/gaps.plot', treatedUnits[i], '.png', sep = ''), width = 10, height = 5)
    
    
    
  }
  
}




