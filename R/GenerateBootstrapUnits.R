#' GenerateBootstrapUnits
#' 
#' A resampling method for geographic quasi-experimental design in which lower geographic units are resampled through bootstrapping
#' so that more control units are created. In synthetic controls, this leads to a better pre-period fit and increased power in the 
#' test statistic.
#' 
#' @param lower_geog_frame - a data frame that contains four columns: 1. lower level geography; 2. higher level geography; 3. date aggregation; 4. count column
#' @param placebo_units - this is a list containing the names of the higher level geography units that did not receive a treatment
#' @param number_drawn - a series of integers that add randomness to the draws of the lower level geography. A number is randomly selected for each new unit created from this list, and this determines the number of draws from the lower level geography.
#' @param number_reps - an integer specifying how many new units should be created for every existing control unit
#' @param lower_geog_colname - a string containing the column name for the lower level geography
#' @param higher_geog_colname - a string containing the column name for the higher level geography
#' @param date_colname - a string containing the column name for the date field
#' @param aggregate_colname - a string containing the column name for the outcome field
#' @param seed - an integer to set the random seed


GenerateBootstrapUnits = function(lower_geog_frame, 
                                    placebo_units, 
                                    number_drawn, 
                                    number_reps,
                                    lower_geog_colname,
                                    higher_geog_colname,
                                    date_colname,
                                    aggregate_colname,
                                    seed = 12345) {
  
  suppressMessages(library(lazyeval))
  suppressMessages(library(dplyr))
  
  set.seed(seed = seed)
  bootstrap_df = data.frame()
  
  for (i in 1:length(placebo_units)) {
    
    placebo_df_old = bootstrap_df
    
    # Admittedly, there is a lot going on here
    sub_unit = sort(unlist(unique(dplyr::select(subset(lower_geog_frame, get(higher_geog_colname) == placebo_units[i]), !!(lower_geog_colname)))))
    new_district_new = data.frame()
    
    for (n in 1:number_reps) {
      
      new_district_old    = new_district_new
      draw                = sample(number_drawn, 1)
      bootstrapped_beats  = sample(sub_unit, draw, replace = T)
      filled_new          = data.frame()
      
      for (b in 1:draw) {
        
        filled_old = filled_new
        filled_tmp = subset(lower_geog_frame, get(lower_geog_colname) == bootstrapped_beats[b])
        filled_new = rbind(filled_old, filled_tmp)
      }
      
      new_district_tmp                  = filled_new %>%
        group_by_(higher_geog_colname, date_colname) %>%
        summarise_(aggregate_colname = (~sum(aggregate_colname)) %>%
                     interp(aggregate_colname = as.name(aggregate_colname))
        )
      names(new_district_tmp) = c(higher_geog_colname, date_colname, aggregate_colname)
      #
      
      
      new_district_tmp[[higher_geog_colname]] = as.character(paste(new_district_tmp[[higher_geog_colname]], n, sep = '_'))
      new_district_new                        = rbind(data.frame(new_district_tmp), new_district_old)
      
    }
    
    
    bootstrap_df = rbind(new_district_new, placebo_df_old)
    
  }
  
  return(bootstrap_df)
  
  
}
