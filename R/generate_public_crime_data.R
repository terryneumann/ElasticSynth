#' Generate public crime data for use with ElasticSynth
#' 
#' @export
generate_public_crime_data <- function() {
  library(dplyr)
  library(RSocrata)
  library(zoo)
  library(lubridate)
  library(tidyr)
  
  fbi_code <- "'05'" # 01A = Murder, 02 = CSA, 03 = Robbery, 04A = Assault, 04B = Battery
  # Pull all data with FBI Code less than 05
  url <- sprintf("https://data.cityofchicago.org/resource/6zsd-86xi.json?$select=*&$where=fbi_code<%s", fbi_code)
  violent_felonies <- read.socrata(url)
  violent_felonies$date_clean <- as.Date(as.POSIXct(violent_felonies$date, format = '%Y-%m-%d %H:%M:%S'))
  violent_felonies$yearmon <- as.Date(as.yearmon(violent_felonies$date_clean))
  
  violent_felonies_by_district <- violent_felonies %>%
    group_by(yearmon, district, .drop=F) %>%
    summarise(countcrimes = length(id)) %>%
    complete(yearmon, district) %>%
    ungroup() %>%
    mutate(fbi_code = 'VF') %>%
    bind_rows(
      violent_felonies %>%
        # define homicide as homicide + second degree homicide
        mutate(fbi_code = ifelse(fbi_code == '01B', '01A', fbi_code)) %>%
        group_by(yearmon, district, fbi_code, .drop =F) %>%
        summarise(countcrimes = length(id)) %>%
        complete(yearmon, district, fbi_code)
    ) %>%
    filter(yearmon >= '2010-03-01' & yearmon <= '2018-12-01') %>%
    # --- period is six month intervals, Period 1 = March through August, Period 2 = September through February. 
    # --- Except for 2017 due to staggered release of SDSCs and introduction of Tier 2 SDSCs before the end of period 2
    mutate(period = ifelse(month(yearmon) %in% 3:8, 1, 2),
           period = ifelse(month(yearmon) %in% 1:2, 
                           as.numeric(paste(year(yearmon) - 1, period, sep = '.')),
                           as.numeric(paste(year(yearmon), period, sep = '.'))),
           period = ifelse(month(yearmon) %in% 1:2 & year(yearmon) == 2018, period + 0.9, period),
           period = ifelse(month(yearmon) == 2 & year(yearmon) == 2017, period + 0.9, period),
           period = as.numeric(as.character(period)),
           period = as.factor(period),
           district = as.factor(district),
           fbi_code = as.factor(fbi_code)) %>%
    group_by(period, district, fbi_code, .drop=F) %>%
    summarise(countcrimes = mean(countcrimes)) %>%
    ungroup() %>%
    # --- generate time variable
    arrange(district, fbi_code, period) %>%
    group_by(district, fbi_code) %>%
    mutate(countcrimes = ifelse(is.nan(countcrimes), 0, countcrimes),
           time = 1:length(period))
  
  treated_units <- subset(violent_felonies_by_district, district %in% c('006', '007', '009', '010', '011', '015'))
  donor_units <- subset(violent_felonies_by_district, !(district %in%  c('006', '007', '009', '010', '011', '015', '031')))
  return(list(treated_units = treated_units,
              donor_units = donor_units))
}
