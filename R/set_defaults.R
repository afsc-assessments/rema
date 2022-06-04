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

set_defaults <- function(input) {

  data = input$data
  par = input$par
  map = input$map

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

  # par$dummy <- 0
  par$log_PE <- rep(1, length(unique(data$pointer_PE_biomass)))
  par$log_q <- rep(1, length(unique(data$pointer_q_cpue)))
  # use linear interpolation to initiate starting values if none are supplied.
  # rule = 2 means the value at the closest data extreme is used
  par$log_biomass_pred <- log(apply(X = data.frame(data$biomass_obs),
                                    MARGIN = 2,
                                    FUN = zoo::na.approx, maxgap = 100, rule = 2))

  # map$log_PE <- fill_vals(par$log_PE, NA)
  # map$log_q <- fill_vals(par$log_q, NA)
  # map$log_biomass_pred <- fill_vals(par$log_biomass_pred, NA)

  input$data = data
  input$par = par
  input$map = map
  return(input)
}
