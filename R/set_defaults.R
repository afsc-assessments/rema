# default REMA data, par, and map list objects for TMB

# default assumptions:
# 1) every biomass stratum gets its own PE. observation error for surveys is lognormal
# 2) if CPUE survey data are included, there are the same number of biomass
# and CPUE survey strata
# 3) every CPUE stratum gets its own q
# 4) no likelihood weights, penalties, or priors
# 5) log_PE and log_q parameters: initial values = 1
# 6) log_pred_biomass random effects: initial values determined using linear
# approximation on biomass_dat observations

set_defaults <- function(input, admb_re = NULL) {

  data = input$data
  par = input$par
  map = input$map
  random = input$random

  data$obs_error_type <- 0

  data$pointer_PE_biomass <- (1:ncol(data$biomass_obs))-1
  data$pointer_biomass_cpue_strata <- (1:ncol(data$biomass_obs))-1
  data$pointer_q_cpue <- (1:ncol(data$cpue_obs))-1
  data$pointer_extra_biomass_cv <- rep(1, ncol(data$biomass_obs))-1
  data$pointer_extra_cpue_cv <- rep(1, ncol(data$cpue_obs))-1

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

  data$extra_biomass_cv <- 0
  data$extra_cpue_cv <- 0

  data$tau_biomass_upper <- rep(1.5, length(unique(data$pointer_extra_biomass_cv)))
  data$tau_cpue_upper <- rep(1.5, length(unique(data$pointer_extra_cpue_cv)))

  # process error and scaling parameters
  par$log_PE <- rep(1, length(unique(data$pointer_PE_biomass)))
  par$log_q <- rep(1, length(unique(data$pointer_q_cpue)))

  # tweedie parameters starting vals
  # tweedie_p_init = 1.6 # based on Dave Warton (ref?)
  # tweedie_p_upper = 2 # gamma
  # tweedie_p_lower = 1 # zero inflated poisson
  # log((tweedie_p_init - tweedie_p_lower) / (tweedie_p_upper - tweedie_p_init))
  # 0.4054651

  if(data$multi_survey == 0) {
    par$logit_tweedie_p <- 0.4054651
  } else {
    par$logit_tweedie_p <- rep(0.4054651, 2)
  }

  # extra CV for biomass or CPUE survey (-Inf = 0 in logit space)
  par$logit_tau_biomass <- rep(-Inf, length(unique(data$pointer_extra_biomass_cv)))
  par$logit_tau_cpue <- rep(-Inf, length(unique(data$pointer_extra_cpue_cv)))

  # if admb_re is provided using read_admb_re(), use log biomass predictions as
  # initial values for the model. if not, use linear interpolation to initiate
  # starting values if none are supplied. rule = 2 means the value at the
  # closest data extreme is used. for the purposes of linear interpolation,
  # remove the very small values artificially generated for zero values

  if(!is.null(admb_re)) {
    par$log_biomass_pred <- admb_re$init_log_biomass_pred
  } else {
    tmp <- data$biomass_obs
    tmp[tmp < 0.01] <- NA # just for starting values
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

  map$logit_tweedie_p <- fill_vals(par$logit_tweedie_p, NA)

  map$logit_tau_biomass <- fill_vals(par$logit_tau_biomass, NA)

  map$logit_tau_cpue <- fill_vals(par$logit_tau_cpue, NA)

  map$log_biomass_pred <- as.factor(1:length(map$log_biomass_pred))

  # by default random = log_pred_biomass
  random = 'log_biomass_pred'

  input$data = data
  input$par = par
  input$map = map
  input$random = random

  return(input)
}
