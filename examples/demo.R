### Generate Crime Data
crime_data <- generate_public_crime_data()
treated_units <- crime_data$treated_units
donor_units <- crime_data$donor_units
#treated_units <- read.csv('~/treated_units.csv')
#donor_units <- read.csv('~/donor_units.csv')
head(treated_units, n = 30)
head(donor_units, n = 30)

#### Run Synthetic Controls
measure_vars = c('01A', '02', '03', '04A', '04B', 'VF')
nvar = length(measure_vars)
pre_list = list(1:14, 1:14, 1:14, 1:14, 1:14, 1:14)
post_list = list(15:18, 15:18, 15:18, 15:18, 15:18, 15:18)
start_cv_list = rep(9, nvar)
end_cv_list = rep(14, nvar)
cv_step_list = rep(2, nvar)
time_unit_list = rep('6 months', nvar)
end_month_list = rep('2018-09-01', nvar)
start_month_list = rep('2010-03-01', nvar)
unit_col = 'district'
time_col = 'time'
measure_col = 'fbi_code'
value_col = 'countcrimes'
lower_limit_weights = 0
upper_limit_weights = Inf
# to-do, clean this up
max_pre_month = as.Date('2016-09-01')
end_month = as.Date('2018-09-01')
start_month = as.Date('2010-03-01')


mod <- ElasticSynth$new(
  treated_units = treated_units,
  donor_units = donor_units,
  pre_list = pre_list,
  post_list = post_list,
  measure_vars = measure_vars,
  measure_col = measure_col,
  unit_col = unit_col,
  time_col = time_col,
  value_col = value_col,
  start_cv_list = start_cv_list,
  end_cv_list = end_cv_list,
  cv_step_list = cv_step_list,
  lower_limit_weights = lower_limit_weights,
  upper_limit_weights = upper_limit_weights,
  start_month_list = start_month_list,
  end_month_list = end_month_list,
  time_unit_list = time_unit_list,
  max_pre_month = max_pre_month,
  start_month = start_month,
  end_month = end_month,
  placebo = T,
  verbose = T)

cv_results <- mod$cv(type = 'individual')
weights <- mod$generate_weights(cv_results = cv_results)
fit <- mod$fit(cv_results = cv_results, weights = weights)
plots <- mod$plot(series_frame = fit$series_frame)

