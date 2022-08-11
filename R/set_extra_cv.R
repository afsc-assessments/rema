# user-defined options for estimating additional cv on biomass or cpue survey
# observations
set_extra_cv <- function(input, extra_biomass_cv, extra_cpue_cv) {

  data = input$data
  par = input$par
  map = input$map

  if((data$multi_survey == 0) & !is.null(extra_cpue_cv)) {
    warning("User has defined inputs for extra_cpue_cv, but the model is in single survey mode. Define multi_survey = 1 in order to estimate additional CPUE observation error.")
  }

  if(!is.null(extra_biomass_cv)) {
    if(!extra_biomass_cv$assumption %in% c('none', 'extra_cv')) {
      stop("User has defined an assumption for extra_biomass_cv other than 'none' or 'extra_cv'.")
    }
    if(extra_biomass_cv$assumption == 'extra_cv') {
      data$extra_biomass_cv <- 1
      map$logit_tau_biomass <- as.factor(1:length(map$logit_tau_biomass))
    }
    if(extra_biomass_cv$assumption == 'none') {
      extra_biomass_cv <- NULL
    }
  }
  if(!is.null(extra_cpue_cv)){
    if(!extra_cpue_cv$assumption %in% c('none', 'extra_cv')) {
      stop("User has defined an assumption for extra_cpue_cv other than 'none' or 'extra_cv'.")
    }
    if(extra_cpue_cv$assumption == 'extra_cv') {
      data$extra_cpue_cv <- 1
      map$logit_tau_cpue <- as.factor(1:length(map$logit_tau_cpue))
    }
    if(extra_cpue_cv$assumption == 'none') {
      extra_cpue_cv <- NULL
    }
  }

  # user defined index for extra biomass CV
  if(!is.null(extra_biomass_cv$pointer_extra_biomass_cv)) {
    if(length(extra_biomass_cv$pointer_extra_biomass_cv) != ncol(data$biomass_obs)) stop("Length of extra_biomass_cv$pointer_extra_biomass_cv must equal the number of biomass survey strata (e.g. length(unique(admb_re$biomass_dat$strata)). Please see extra_biomass_cv details in ?prepare_rema_input.")
    if(any(!is.numeric(extra_biomass_cv$pointer_extra_biomass_cv))) stop("extra_biomass_cv$pointer_extra_biomass_cv must be a vector of integer values starting at 1 with a vector length equal the number of biomass survey strata (e.g. length(unique(admb_re$biomass_dat$strata)). Please see extra_biomass_cv details in ?prepare_rema_input.")
    extra_biomass_cv$pointer_extra_biomass_cv <- as.integer(extra_biomass_cv$pointer_extra_biomass_cv)
    if(any(!is.integer(extra_biomass_cv$pointer_extra_biomass_cv))) stop("extra_biomass_cv$pointer_extra_biomass_cv must be a vector of integer values starting at 1 with a vector length equal the number of biomass survey strata (e.g. length(unique(admb_re$biomass_dat$strata)). Please see extra_biomass_cv details in ?prepare_rema_input.")
    data$extra_biomass_cv <- 1
    data$pointer_extra_biomass_cv <- (extra_biomass_cv$pointer_extra_biomass_cv)-1 # TMB started indexing at 0
    data$tau_biomass_upper <- rep(1.5, length(unique(data$pointer_extra_biomass_cv)))
    par$logit_tau_biomass <- rep(1, length(unique(data$pointer_extra_biomass_cv)))
    map$logit_tau_biomass <- par$logit_tau_biomass
    map$logit_tau_biomass <- as.factor(1:length(map$logit_tau_biomass))
  }

  # user defined index for extra cpue CV
  if(!is.null(extra_cpue_cv$pointer_extra_cpue_cv)) {
    if(length(extra_cpue_cv$pointer_extra_cpue_cv) != ncol(data$cpue_obs)) stop("Length of extra_cpue_cv$pointer_extra_cpue_cv must equal the number of cpue survey strata (e.g. length(unique(admb_re$cpue_dat$strata)). Please see extra_cpue_cv details in ?prepare_rema_input.")
    if(any(!is.numeric(extra_cpue_cv$pointer_extra_cpue_cv))) stop("extra_cpue_cv$pointer_extra_cpue_cv must be a vector of integer values starting at 1 with a vector length equal the number of cpue survey strata (e.g. length(unique(admb_re$cpue_dat$strata)). Please see extra_cpue_cv details in ?prepare_rema_input.")
    extra_cpue_cv$pointer_extra_cpue_cv <- as.integer(extra_cpue_cv$pointer_extra_cpue_cv)
    if(any(!is.integer(extra_cpue_cv$pointer_extra_cpue_cv))) stop("extra_cpue_cv$pointer_extra_cpue_cv must be a vector of integer values starting at 1 with a vector length equal the number of cpue survey strata (e.g. length(unique(admb_re$cpue_dat$strata)). Please see extra_cpue_cv details in ?prepare_rema_input.")
    data$extra_cpue_cv <- 1
    data$pointer_extra_cpue_cv <- (extra_cpue_cv$pointer_extra_cpue_cv)-1 # TMB started indexing at 0
    data$tau_cpue_upper <- rep(1.5, length(unique(data$pointer_extra_cpue_cv)))
    par$logit_tau_cpue <- rep(1, length(unique(data$pointer_extra_cpue_cv)))
    map$logit_tau_cpue <- par$logit_tau_cpue
    map$logit_tau_cpue <- as.factor(1:length(map$logit_tau_cpue))
  }

  # user defined initial values for logit_tau_biomass parameters
  if(!is.null(extra_biomass_cv$initial_pars)) {
    if(length(extra_biomass_cv$initial_pars) != length(par$logit_tau_biomass)) stop("extra_biomass_cv$initial_pars must be a vector with length equal to the number of extra biomass cv (logit_tau_biomass) parameters in the model.  Input values must be numeric and users should specify logit_tau_biomass initial values in logit space.")
    if(any(!is.numeric(extra_biomass_cv$initial_pars))) stop("extra_biomass_cv$initial_pars must be a vector with length equal to the number of extra biomass cv (logit_tau_biomass) parameters in the model.  Input values must be numeric and users should specify logit_tau_biomass initial values in logit space.")
    par$logit_tau_biomass <- extra_biomass_cv$initial_pars
  }

  # user defined initial values for logit_tau_cpue parameters
  if(!is.null(extra_cpue_cv$initial_pars)) {
    if(length(extra_cpue_cv$initial_pars) != length(par$logit_tau_cpue)) stop("extra_cpue_cv$initial_pars must be a vector with length equal to the number of extra cpue cv (logit_tau_cpue) parameters in the model.  Input values must be numeric and users should specify logit_tau_cpue initial values in logit space.")
    if(any(!is.numeric(extra_cpue_cv$initial_pars))) stop("extra_cpue_cv$initial_pars must be a vector with length equal to the number of extra cpue cv (logit_tau_cpue) parameters in the model.  Input values must be numeric and users should specify logit_tau_cpue initial values in logit space.")
    par$logit_tau_cpue <- extra_cpue_cv$initial_pars
  }

  # user defined fixed values of logit_tau_biomass
  if(!is.null(extra_biomass_cv$fix_pars)) {
    if(any(!is.numeric(extra_biomass_cv$fix_pars))) stop("extra_biomass_cv$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tau_biomass to be fixed.")
    extra_biomass_cv$fix_pars <- as.integer(extra_biomass_cv$fix_pars)
    if(any(!is.integer(extra_biomass_cv$fix_pars))) stop("extra_biomass_cv$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tau_biomass to be fixed.")

    # create logit_tau_biomass map for fixed/estimated parameters
    tmp_logit_tau_biomass <- par$logit_tau_biomass
    icounter <- 1
    for(i in 1:length(tmp_logit_tau_biomass)) {
      if(i %in% extra_biomass_cv$fix_pars) {
        tmp_logit_tau_biomass[i] <- NA
      } else {
        tmp_logit_tau_biomass[i] <- icounter
        icounter <- icounter + 1
      }
    }
    map$logit_tau_biomass <- as.factor(tmp_logit_tau_biomass)
  }

  # user defined fixed values of logit_tau_cpue
  if(!is.null(extra_cpue_cv$fix_pars)) {
    if(any(!is.numeric(extra_cpue_cv$fix_pars))) stop("extra_cpue_cv$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tau_cpue to be fixed.")
    extra_cpue_cv$fix_pars <- as.integer(extra_cpue_cv$fix_pars)
    if(any(!is.integer(extra_cpue_cv$fix_pars))) stop("extra_cpue_cv$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tau_cpue to be fixed.")

    # create logit_tau_cpue map for fixed/estimated parameters
    tmp_logit_tau_cpue <- par$logit_tau_cpue
    icounter <- 1
    for(i in 1:length(tmp_logit_tau_cpue)) {
      if(i %in% extra_cpue_cv$fix_pars) {
        tmp_logit_tau_cpue[i] <- NA
      } else {
        tmp_logit_tau_cpue[i] <- icounter
        icounter <- icounter + 1
      }
    }
    map$logit_tau_cpue <- as.factor(tmp_logit_tau_cpue)
  }

  input$data = data
  input$par = par
  input$map = map

  return(input)
}
