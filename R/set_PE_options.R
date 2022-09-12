# user-defined process error options, either to default values or settings specified by
# user in the PE_options list
set_PE_options <- function(input, PE_options) {

  data = input$data
  par = input$par
  map = input$map

  # user defined index for PE estimation (e.g. there are 3 strata but user
  # only wants to estimate 1 PE, pointer_PE_biomass = c(1, 1, 1))
  if(!is.null(PE_options$pointer_PE_biomass)) {
    if(length(PE_options$pointer_PE_biomass) != ncol(data$biomass_obs)) stop("Length of PE_options$pointer_PE_biomass must equal the number of biomass survey strata (e.g. length(unique(admb_re$biomass_dat$strata)). Please see PE_options details in ?prepare_rema_input.")
    PE_options$pointer_PE_biomass <- as.integer(PE_options$pointer_PE_biomass)
    if(!any(is.integer(PE_options$pointer_PE_biomass))) stop("PE_options$pointer_PE_biomass must be a vector of integer values starting at 1 with a vector length equal the number of biomass survey strata (e.g. length(unique(admb_re$biomass_dat$strata)). Please see PE_options details in ?prepare_rema_input.")
    data$pointer_PE_biomass <- (PE_options$pointer_PE_biomass)-1 # TMB started indexing at 0
    par$log_PE <- rep(1, length(unique(data$pointer_PE_biomass)))
    map$log_PE <- par$log_PE
    map$log_PE <- as.factor(1:length(map$log_PE))
  }

  # user defined initial values for log_PE parameters
  if(!is.null(PE_options$initial_pars)) {
    if(length(PE_options$initial_pars) != length(par$log_PE)) stop("PE_options$initial_pars must be a vector with length equal to the number of PE parameters in the model.  Input values must be numeric and users should specify log_PE initial values in log space.")
    if(!any(is.numeric(PE_options$initial_pars))) stop("PE_options$initial_pars must be a vector with length equal to the number of PE parameters in the model. Input values must be numeric and users should specify log_PE initial values in log space.")
    par$log_PE <- PE_options$initial_pars
  }

  # user defined fixed values of log_PE. note that this is not recommended.
  if(!is.null(PE_options$fix_pars)) {
    warning("Are you sure you want to fix process error parameters? This is not recommended but may be useful for sensitivity analysis or troubleshooting. User can fix log_PE value to something other than the default of log_PE = 1 using the PE_options$initial_pars setting.")

    PE_options$fix_pars <- as.integer(PE_options$fix_pars)
    if(!any(is.integer(PE_options$fix_pars))) stop("PE_options$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the log_PE to be fixed. Please see PE_options details in ?prepare_rema_input.")

    # create log_PE map for fixed/estimated parameters
    tmp_log_PE <- par$log_PE
    icounter <- 1
    for(i in 1:length(tmp_log_PE)) {
      if(i %in% PE_options$fix_pars) {
        tmp_log_PE[i] <- NA
      } else {
        tmp_log_PE[i] <- icounter
        icounter <- icounter + 1
      }
    }
    map$log_PE <- as.factor(tmp_log_PE)
  }

  # PE penalties and priors
  if(!is.null(PE_options$penalty_options)) {

    if(!PE_options$penalty_options %in% c('none', 'wt', 'squared_penalty', 'normal_prior')) {
      stop("User has defined an invalid input for PE_options$penalty_options. Please check the Details section in ?prepare_rema_input.")
    }

    if(PE_options$penalty_options == 'wt') {

      if(length(PE_options$penalty_values) > 1) stop("User has selected PE_options$penalty_options = 'wt' but has supplied more than one value in PE_options$penalty_values. Please provide a single, numeric value >= 0.")
      if(!any(is.numeric(PE_options$penalty_values))) stop("User has selected PE_options$penalty_options = 'wt' but has supplied a non-numeric value in PE_options$penalty_values. Please provide a single, numeric value >= 0.")

      data$PE_penalty_type <- 1
      data$wt_PE <- PE_options$penalty_values
    }

    if(PE_options$penalty_options == 'squared_penalty') {

      if(length(PE_options$penalty_values) != length(par$log_PE)) stop("User has selected PE_options$penalty_options = 'squared_penalty' but has supplied an incorrect number of values to PE_options$penalty_values. Please provide a vector of numeric values with length equal to the number of PE parameters.")
      if(!any(is.numeric(PE_options$penalty_values))) stop("User has selected PE_options$penalty_options = 'squared_penalty' but has supplied non-numeric values to PE_options$penalty_values. Please provide a vector of numeric values with length equal to the number of PE parameters.")

      data$PE_penalty_type <- 2
      data$squared_penalty_PE <- PE_options$penalty_values
    }

    if(PE_options$penalty_options == 'normal_prior') {

      if(length(PE_options$penalty_values) != 2 * length(par$log_PE)) stop("User has selected PE_options$penalty_options = 'normal_prior' but has supplied an incorrect number of values to PE_options$penalty_values. See PE_options details in ?prepare_rema_input for examples of valid inputs.")
      if(!any(is.numeric(PE_options$penalty_values))) stop("User has selected PE_options$penalty_options = 'normal_prior' but has supplied non-numeric values to PE_options$penalty_values. See PE_options details in ?prepare_rema_input for examples of valid inputs.")

      data$PE_penalty_type <- 3
      data$pmu_log_PE <- PE_options$penalty_values[seq(1,length(PE_options$penalty_values),2)]
      data$psig_log_PE <- PE_options$penalty_values[seq(0,length(PE_options$penalty_values),2)]
    }
  }
  input$data = data
  input$par = par
  input$map = map

  return(input)
}
