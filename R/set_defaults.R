# default REMA data, par, and map list objects for TMB

# default assumptions:
# 1) every biomass stratum gets its own PE
# 2) if CPUE survey data are included, there are the same number of biomass
# and CPUE survey strata
# 3) every CPUE stratum gets its own q
# 4) no likelihood weights, penalties, or priors
# 5) log_PE and log_q parameters: initial values = 1
# 6) log_pred_biomass random effects: initial values determined using linear
# approximation on biomass_dat observations

set_defaults <- function(input, re_dat = NULL) {

  data = input$data
  par = input$par
  map = input$map
  random = input$random

  data$pointer_PE_biomass <- (1:ncol(data$biomass_obs))-1
  data$pointer_q_biomass <- (1:ncol(data$biomass_obs))-1
  data$pointer_q_cpue <- (1:ncol(data$cpue_obs))-1

  data$wt_biomass <- 1
  data$wt_cpue <- 1

  data$PE_penalty_type <- 0
  data$wt_PE <- 1
  data$squared_penalty_PE <- NA
  data$pmu_log_PE <- NA
  data$psig_log_PE <- NA

  data$q_penalty_type <- 0
  data$pmu_log_q <- NA
  data$psig_log_q <- NA

  par$log_PE <- rep(1, length(unique(data$pointer_PE_biomass)))
  par$log_q <- rep(1, length(unique(data$pointer_q_cpue)))

  # if re_dat is provided, use log biomass predictions as initial values for the
  # model. if not, use linear interpolation to initiate starting values if none
  # are supplied. rule = 2 means the value at the closest data extreme is used.
  # for the purposes of linear interpolation, remove the very small values
  # artificially generated for zero values

  if(!is.null(re_dat)) {
    par$log_biomass_pred <- re_dat$init_log_biomass_pred
  } else {
  tmp <- data$biomass_obs
  tmp[tmp < 0.0001] <- NA
  par$log_biomass_pred <- log(apply(X = tmp,
                                    MARGIN = 2,
                                    FUN = zoo::na.approx, maxgap = 100, rule = 2))
  }

  # set map
  map <- par
  map$log_PE <- as.factor(1:length(map$log_PE))

  if(data$multi_survey == 0) {
    map$log_q <- fill_vals(par$log_q, NA)
  } else {
    map$log_q <- as.factor(1:length(map$log_q))
  }

  # FLAG I don't think the parameter factor order level matters here but if it
  # does do I need a byrow = TRUE
  map$log_biomass_pred <- as.factor(1:length(map$log_biomass_pred))

  # by default random = log_pred_biomass
  random = 'log_biomass_pred'

  input$data = data
  input$par = par
  input$map = map
  input$random = random

  return(input)
}
