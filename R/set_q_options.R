# user-defined scaling parameters options. settings specified by user in the
# q_options list
set_q_options <- function(input, q_options) {

  data = input$data
  par = input$par
  map = input$map

  if((data$multi_survey == 1) &
     ncol(input$data$biomass_obs) != ncol(input$data$cpue_obs) &
     is.null(q_options$pointer_biomass_cpue_strata)) {
    stop("Strata definitions differ for the biomass and CPUE surveys. The user must define q_options$pointer_biomass_cpue_strata, a vector with length = number of biomass survey strata and unique values = the number of CPUE survey strata. See q_options details in ?prepare_rema_input for examples of valid inputs.")
  }

  # user defined index for q estimation by CPUE survey strata (e.g. there are
  # 2 strata but user only wants to estimate 1 q, pointer_q_cpue = c(1, 1)).
  # FLAG - I can't actually think of a user scenario when this would be
  # desirable but thought I would build it in now instead of later
  if(!is.null(q_options$pointer_q_cpue)) {
    if(length(q_options$pointer_q_cpue) != ncol(data$cpue_obs)) stop("Length of q_options$pointer_q_cpue must equal the number of CPUE survey strata (e.g. length(unique(admb_re$cpue_dat$strata)). Please see q_options details in ?prepare_rema_input.")
    if(any(!is.numeric(q_options$pointer_q_cpue))) stop("q_options$pointer_q_cpue must be a vector of integer values starting at 1 with a vector length equal the number of CPUE survey strata (e.g. length(unique(admb_re$cpue_dat$strata)). Please see q_options details in ?prepare_rema_input.")
    q_options$pointer_q_cpue <- as.integer(q_options$pointer_q_cpue)
    if(any(!is.integer(q_options$pointer_q_cpue))) stop("q_options$pointer_q_cpue must be a vector of integer values starting at 1 with a vector length equal the number of CPUE survey strata (e.g. length(unique(admb_re$cpue_dat$strata)). Please see q_options details in ?prepare_rema_input.")
    data$pointer_q_cpue <- (q_options$pointer_q_cpue)-1 # TMB started indexing at 0
    par$log_q <- rep(1, length(unique(data$pointer_q_cpue)))
    map$log_q <- par$log_q
    map$log_q <- as.factor(1:length(map$log_q))
    # FLAG ASK PH:
    # if(length(unique(q_options$pointer_q_cpue)) != ncol(data$cpue_obs)) warning(paste0("The recommended model configuration is to estimate one scaling parameter (log_q) for each CPUE survey stratum. Your inputs for q_options$pointer_q_cpue currently specify ", length(unique(q_options$pointer_q_cpue)), " log_q for ", ncol(data$cpue_obs), " CPUE survey strata."))
  }

  # user defined index for biomass predictions by CPUE strata (e.g. there are 3
  # biomass strata but only two CPUE strata. if pointer_biomass_cpue_strata =
  # c(1, 1, 2)), this means the first two biomass strata correspond to the first
  # CPUE survey stratum, and the third biomass strata corresponds to second CPUE
  # survey stratum.
  if(!is.null(q_options$pointer_biomass_cpue_strata)) {
    if(length(q_options$pointer_biomass_cpue_strata) != ncol(data$biomass_obs)) stop("Length of q_options$pointer_biomass_cpue_strata must equal the number of biomass survey strata (e.g., length(unique(admb_re$biomass_dat$strata)). Please see q_options details in ?prepare_rema_input.")
    if(length(unique(q_options$pointer_biomass_cpue_strata)) != ncol(data$cpue_obs)) stop("Length of unique values in q_options$pointer_biomass_cpue_strata must equal the number of CPUE survey strata (e.g., length(unique(admb_re$cpue_dat$strata)). Please see q_options details in ?prepare_rema_input.")
    if(any(!is.numeric(q_options$pointer_biomass_cpue_strata))) stop("q_options$pointer_biomass_cpue_strata must be a vector of integer values starting at 1 with a vector length equal the number of biomass survey strata (e.g., length(unique(admb_re$biomass_dat$strata)). Please see q_options details in ?prepare_rema_input.")
    q_options$pointer_biomass_cpue_strata <- as.integer(q_options$pointer_biomass_cpue_strata)
    if(any(!is.integer(q_options$pointer_biomass_cpue_strata))) stop("q_options$pointer_biomass_cpue_strata must be a vector of integer values starting at 1 with a vector length equal the number of biomass survey strata (e.g., length(unique(admb_re$biomass_dat$strata)). Please see q_options details in ?prepare_rema_input.")
    data$pointer_biomass_cpue_strata <- (q_options$pointer_biomass_cpue_strata)-1 # TMB started indexing at 0
  }

  # user defined initial values for log_q parameters
  if(!is.null(q_options$initial_pars)) {
    if(length(q_options$initial_pars) != length(par$log_q)) stop("q_options$initial_pars must be a vector with length equal to the number of q parameters in the model.  Input values must be numeric and users should specify log_q initial values in log space.")
    if(any(!is.numeric(q_options$initial_pars))) stop("q_options$initial_pars must be a vector with length equal to the number of q parameters in the model. Input values must be numeric and users should specify log_q initial values in log space.")
    par$log_q <- q_options$initial_pars
  }

  # user defined fixed values of log_q.
  if(!is.null(q_options$fix_pars)) {
    if(any(!is.numeric(q_options$fix_pars))) stop("q_options$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the log_q to be fixed.")
    q_options$fix_pars <- as.integer(q_options$fix_pars)
    if(any(!is.integer(q_options$fix_pars))) stop("q_options$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the log_q to be fixed.")

    # create log_q map for fixed/estimated parameters
    tmp_log_q <- par$log_q
    icounter <- 1
    for(i in 1:length(tmp_log_q)) {
      if(i %in% q_options$fix_pars) {
        tmp_log_q[i] <- NA
      } else {
        tmp_log_q[i] <- icounter
        icounter <- icounter + 1
      }
    }
    map$log_q <- as.factor(tmp_log_q)
  }

  # q penalties and priors
  if(!is.null(q_options$penalty_options)) {

    if(!q_options$penalty_options %in% c('none', 'normal_prior')) {
      stop("User has defined an invalid input for q_options$penalty_options. Please check q_options details in ?prepare_rema_input.")
    }

    if(q_options$penalty_options == 'normal_prior') {

      if(length(q_options$penalty_values) != 2 * length(par$log_q)) stop("User has selected q_options$penalty_options = 'normal_prior' but has supplied an incorrect number of values to q_options$penalty_values. See q_options details in ?prepare_rema_input for examples of valid inputs.")
      if(any(!is.numeric(q_options$penalty_values))) stop("User has selected q_options$penalty_options = 'normal_prior' but has supplied non-numeric values to q_options$penalty_values. See q_options details in ?prepare_rema_input for examples of valid inputs.")

      data$q_penalty_type <- 1
      data$pmu_log_q <- q_options$penalty_values[seq(1,length(q_options$penalty_values),2)]
      data$psig_log_q <- q_options$penalty_values[seq(0,length(q_options$penalty_values),2)]
    }
  }

  input$data = data
  input$par = par
  input$map = map

  return(input)
}
