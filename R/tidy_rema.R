#' Tidy REMA model output
#'
#' Takes outputs from \code{\link{fit_rema}}, and returns a named list of tidied
#' data.frames that include parameter estimates and standard errors, and derived
#' variables from the model. For more information on "tidy" data, please see
#' \href{https://vita.had.co.nz/papers/tidy-data.pdf}{Wickham 2014}. Some code
#' modified from
#' \href{https://github.com/timjmiller/wham/blob/master/R/par_tables_fn.R}{\code{wham::par_tables_fun}}.
#'
#' @param rema_model list out output from \code{\link{fit_rema}}, which includes
#'   model results but also inputs
#' @param save (optional) logical (T/F) save output data.frames as csvs in
#'   \code{path}. Default = FALSE. NOT YET IMPLEMENTED.
#' @param path (optional) directory path to location where csvs are to be saved
#'   if \code{save = TRUE}. NOT YET IMPLEMENTED.
#' @param alpha_ci (optional) the significance level for generating confidence
#'   intervals. Default = 0.05
#'
#' @return a list with the following items:
#'   \describe{
#'     \item{\code{$parameter_estimates}}{A data.frame of fixed effects
#'     parameters in REMA (e.g. log_PE and log_q) with standard errors and
#'     confidence intervals that have been transformed from log space to natural
#'     space for ease of interpretation. }
#'     \item{\code{$biomass_by_strata}}{A tidy, long format data.frame of model
#'     predicted and observed biomass by biomass survey strata.}
#'     \item{\code{$cpue_by_strata}}{A tidy, long format data.frame of model
#'     predicted and observed CPUE by CPUE survey strata. If REMA is not run in
#'     multi-survey mode, or if CPUE data are not provided, an explanatory
#'     character string with instructions for fitting to CPUE data is returned.}
#'     \item{\code{$biomass_by_cpue_strata}}{A tidy, long format data.frame of
#'     model predicted biomass by CPUE survey strata. Note that observed/summed
#'     biomass observations are not returned in case there are missing values in
#'     one stratum but not another within a given year. This output is reserved
#'     for instances when the number of biomass strata exceeds that of CPUE
#'     survey strata, but the user wants to visualize predicted biomass at the
#'     same resolution as the CPUE predictions. In other scenarios, a character
#'     string is returned explaining the special use case for this object.}
#'     \item{\code{$total_predicted_biomass}}{A tidy, long format data.frame of
#'     total model predicted biomass summed across all biomass survey strata. If
#'     only one stratum is used (i.e. the univariate RE), the predicted values
#'     will be the same as \code{output$biomass_by_strata}.}
#'     \item{\code{$total_predicted_cpue}}{A tidy, long format data.frame of
#'     total model predicted CPUE summed across all CPUE survey strata. If only
#'     one stratum is used (i.e. the univariate RE), the predicted values will
#'     be the same as \code{output$cpue_by_strata}. If The CPUE survey index
#'     provided was defined as not summable in prepare_rema_input(), an
#'     character string will be returned explaining how to change this using the
#'     'sum_cpue_index' in \code{?prepare_rema_input} if appropriate.}
#'   }
#' @export
#' @seealso \code{\link{fit_rema}}
#' @examples
#' \dontrun{
#' # placeholder for example
#' }
tidy_rema <- function(rema_model,
                      save = FALSE, # NOT IMPLEMENTED
                      path = NULL, # NOT IMPLEMENTED
                      alpha_ci = 0.05) {

  # rema_model = m1
  # alpha_ci = 0.05

  if(isFALSE(rema_model$is_sdrep)) {
    stop("Please run fit_rema() with 'do.sdrep = TRUE' in order to get tidied output of estimated and derived variables with standard errors and confidence intervals. See ?fit_rema for details.")
  }

  data <- rema_model$input$data
  sdrep <- rema_model$sdrep

  # parameter table

  # function to get confidence intervals
  get_ci <- function(par, se, alpha = alpha_ci, lo = 0, hi = 1, type = "I", k = 1){

    p = 1 - alpha/2
    ci = par + c(-1, 1) * qnorm(p) * se

    if(type == "I") {
      return(c(se, ci))
    }
    if(type == "exp") {
      return(c(exp(par) * se, exp(ci))) # approximation in natural space is the cv
    }
    if(type == "expit") { # Delta-method: V(lo + (hi-lo)/(1 + exp(-x))) ~ ((hi-lo) * p * (1-p))^2 * V(x)
      p = 1/(1 + exp(- k * par))
      dm.se = k * abs(hi - lo) * p * (1 - p) * se
      return(c(dm.se, lo + (hi-lo)/(1+ exp(-ci))))
    }
  }

  pars <- as.list(sdrep, "Est")
  sd <- as.list(sdrep, "Std")

  # process error parameters
  pe_pars <- NULL

  if(length(pars$log_PE) == 0) {
    warning("No process error parameters were estimated for this model. Output is likely invalid.")
  } else {

    pe_pars <- data.frame(model_name = rema_model$input$model_name,
                          parameter = 'process_error',
               estimate = exp(pars$log_PE))

    pe_ci <- matrix(nrow = 0, ncol = 3)
    for(i in 1:length(pars$log_PE)) {
      pe_ci <- rbind(pe_ci, get_ci(pars$log_PE[i], sd$log_PE[i], type = 'exp'))
    }

    pe_ci <- as.data.frame(pe_ci)
    names(pe_ci) <- c('std_err', 'lci', 'uci')

    pe_pars <- pe_pars %>%
      dplyr::bind_cols(pe_ci)
  }

  # scaling parameter(s) when available
  q_pars <- NULL
  if(data$multi_survey == 1 & length(pars$log_q) > 0) {
    q_pars <- data.frame(model_name = rema_model$input$model_name,
                         parameter = 'scaling_parameter_q',
                         estimate = exp(pars$log_q))

    q_ci <- matrix(nrow = 0, ncol = 3)
    for(i in 1:length(pars$log_q)) {
      q_ci <- rbind(q_ci, get_ci(pars$log_q[i], sd$log_q[i], type = 'exp'))
    }

    q_ci <- as.data.frame(q_ci)
    names(q_ci) <- c('std_err', 'lci', 'uci')

    q_pars <- q_pars %>%
      dplyr::bind_cols(q_ci)
  }

  # power parameter(s) when available
  p_pars <- NULL
  if(any(grepl('logit_tweedie_p', names(pars)))) {
    p_pars <- data.frame(model_name = rema_model$input$model_name,
                         parameter = 'tweedie_p',
                         estimate = 1 + (1 / (1 + exp(-pars$logit_tweedie_p))))


    p_ci <- matrix(nrow = 0, ncol = 3)
    for(i in 1:length(pars$logit_tweedie_p)) {
      p_ci <- rbind(p_ci, get_ci(pars$logit_tweedie_p[i], sd$logit_tweedie_p[i], type = 'expit', lo = 1, hi = 2))
    }

    p_ci <- as.data.frame(p_ci)
    names(p_ci) <- c('std_err', 'lci', 'uci')

    p_pars <- p_pars %>%
      dplyr::bind_cols(p_ci)
  }

  # biomass tau (extra CV) parameter(s) when available
  tau_biomass_pars <- NULL
  if(any(grepl('logit_tau_biomass', names(pars)))) {
    tau_biomass_pars <- data.frame(model_name = rema_model$input$model_name,
                         parameter = 'extra_biomass_cv',
                         estimate = 1.5 / (1 + exp(-pars$logit_tau_biomass)))


    tau_biomass_ci <- matrix(nrow = 0, ncol = 3)
    for(i in 1:length(pars$logit_tau_biomass)) {
      tau_biomass_ci <- rbind(tau_biomass_ci, get_ci(pars$logit_tau_biomass[i], sd$logit_tau_biomass[i], type = 'expit', lo = 0, hi = data$tau_biomass_upper[i]))
    }

    tau_biomass_ci <- as.data.frame(tau_biomass_ci)
    names(tau_biomass_ci) <- c('std_err', 'lci', 'uci')

    tau_biomass_pars <- tau_biomass_pars %>%
      dplyr::bind_cols(tau_biomass_ci)
  }

  # cpue tau (extra CV) parameter(s) when available
  tau_cpue_pars <- NULL
  if(any(grepl('logit_tau_cpue', names(pars)))) {
    tau_cpue_pars <- data.frame(model_name = rema_model$input$model_name,
                                   parameter = 'extra_cpue_cv',
                                   estimate = 1.5 / (1 + exp(-pars$logit_tau_cpue)))


    tau_cpue_ci <- matrix(nrow = 0, ncol = 3)
    for(i in 1:length(pars$logit_tau_cpue)) {
      tau_cpue_ci <- rbind(tau_cpue_ci, get_ci(pars$logit_tau_cpue[i], sd$logit_tau_cpue[i], type = 'expit', lo = 0, hi = data$tau_cpue_upper[i]))
    }

    tau_cpue_ci <- as.data.frame(tau_cpue_ci)
    names(tau_cpue_ci) <- c('std_err', 'lci', 'uci')

    tau_cpue_pars <- tau_cpue_pars %>%
      dplyr::bind_cols(tau_cpue_ci)
  }
  parameter_estimates <- NULL

  if(!is.null(q_pars)) {
    parameter_estimates <- dplyr::bind_rows(pe_pars, q_pars)
  } else {
    parameter_estimates <- pe_pars
  }

  if(!is.null(p_pars) & data$obs_error_type == 1) {
    parameter_estimates <- dplyr::bind_rows(parameter_estimates, p_pars)
  }

  if(!is.null(tau_biomass_pars) & data$extra_biomass_cv == 1) {
    parameter_estimates <- dplyr::bind_rows(parameter_estimates, tau_biomass_pars)
  }

  if(!is.null(tau_cpue_pars) & data$extra_cpue_cv == 1) {
    parameter_estimates <- dplyr::bind_rows(parameter_estimates, tau_cpue_pars)
  }

  # Model predictions:

  # Model estimates of biomass by strata
  ts_biomass_strata <- NULL
  ts_biomass_strata <- tidyr::expand_grid(model_name = rema_model$input$model_name,
                                          strata = colnames(data$biomass_obs),
                                          variable = 'biomass_pred',
                                          year = data$model_yrs) %>%
    dplyr::mutate(log_pred = sdrep$value[names(sdrep$value) == 'log_biomass_pred'],
                  sd_log_pred = sdrep$sd[which(names(sdrep$value) == 'log_biomass_pred')],
                  pred = exp(log_pred),
                  pred_lci = exp(log_pred - qnorm(1 - alpha_ci/2) * sd_log_pred),
                  pred_uci = exp(log_pred + qnorm(1 - alpha_ci/2) * sd_log_pred))

  if(is.null(ts_biomass_strata)) {
    stop("Something went wrong... Please review the ?prepare_rema_input, ?fit_rema, and make sure REMA converged by running check_convergence(my_rema_model).")
  }

  biomass_summary <- NULL

  biomass_summary <- ts_biomass_strata %>%
    dplyr::filter(variable == 'biomass_pred') %>%
    dplyr::left_join(rema_model$input$biomass_dat %>%
                       dplyr::rename(obs = biomass, obs_cv = cv) %>%
                       dplyr::mutate(log_obs = ifelse(obs > 0, log(obs), NA),
                                     sd_log_obs = ifelse(obs > 0, sqrt(log(obs_cv^2 + 1)), NA),
                                     obs_lci = exp(log_obs - qnorm(1 - alpha_ci/2) * sd_log_obs),
                                     obs_uci = exp(log_obs + qnorm(1 - alpha_ci/2) * sd_log_obs)))

  # Model estimates of cpue by strata when available
  ts_cpue_strata <- NULL
  cpue_summary <- NULL
  biomass_by_cpue_strata <- NULL
  total_predicted_cpue <- NULL

  if(data$multi_survey == 0) {
    total_predicted_cpue <- "REMA was fit only to biomass survey data, therefore no CPUE predictions available. If the user has a CPUE survey index and wants to fit to it, please see ?prepare_rema_input() for details."
  }

  if(data$multi_survey == 1) {
    ts_cpue_strata <- tidyr::expand_grid(model_name = rema_model$input$model_name,
                                         strata = colnames(data$cpue_obs),
                                         variable = 'cpue_pred',
                                         year = data$model_yrs) %>%
      dplyr::mutate(log_pred = sdrep$value[names(sdrep$value) == 'log_cpue_pred'],
                    sd_log_pred = sdrep$sd[which(names(sdrep$value) == 'log_cpue_pred')],
                    pred = exp(log_pred),
                    pred_lci = exp(log_pred - qnorm(1 - alpha_ci/2) * sd_log_pred),
                    pred_uci = exp(log_pred + qnorm(1 - alpha_ci/2) * sd_log_pred))

    cpue_summary <- ts_cpue_strata %>%
      dplyr::filter(variable == 'cpue_pred') %>%
      dplyr::left_join(rema_model$input$cpue_dat %>%
                         dplyr::rename(obs = cpue, obs_cv = cv) %>%
                         dplyr::mutate(log_obs = ifelse(obs > 0, log(obs), NA),
                                       sd_log_obs = ifelse(obs > 0, sqrt(log(obs_cv^2 + 1)), NA),
                                       obs_lci = exp(log_obs - qnorm(1 - alpha_ci/2) * sd_log_obs),
                                       obs_uci = exp(log_obs + qnorm(1 - alpha_ci/2) * sd_log_obs)))

    # get total predicted cpue
    if(ncol(data$cpue_obs) == 1) {

      total_predicted_cpue <- ts_cpue_strata %>%
        dplyr::select(model_name, variable, year, pred, pred_lci, pred_uci) %>%
        dplyr::mutate(variable = 'tot_cpue_pred')


       } else if(data$sum_cpue_index == 0) {
        total_predicted_cpue <- "The CPUE survey index provided was defined as not summable in prepare_rema_input(). If the CPUE index is summable (e.g. Relative Population Numbers), please select sum_cpue_index = TRUE in prepare_rema_input(). See ?prepare_rema_input() for more details."

       } else {

         total_predicted_cpue <- tidyr::expand_grid(model_name = rema_model$input$model_name,
                                                    variable = 'tot_cpue_pred',
                                                    year = data$model_yrs) %>%
           dplyr::mutate(log_pred = sdrep$value[names(sdrep$value) == 'log_tot_cpue_pred'],
                         sd_log_pred = sdrep$sd[which(names(sdrep$value) == 'log_tot_cpue_pred')],
                         pred = exp(log_pred),
                         pred_lci = exp(log_pred - qnorm(1 - alpha_ci/2) * sd_log_pred),
                         pred_uci = exp(log_pred + qnorm(1 - alpha_ci/2) * sd_log_pred)) %>%
           dplyr::select(model_name, variable, year, pred, pred_lci, pred_uci)
       }

      # if there are more biomass strata than cpue strata, get
      # predicted biomass by cpue strata for comparison at the same strata level
      if(any(grepl('log_biomass_pred_cpue_strata', unique(names(sdrep$value))))) {

      biomass_by_cpue_strata <- tidyr::expand_grid(model_name = rema_model$input$model_name,
                                           strata = colnames(data$cpue_obs),
                                           variable = 'biomass_pred_cpue_strata',
                                           year = data$model_yrs) %>%
        dplyr::mutate(log_pred = sdrep$value[names(sdrep$value) == 'log_biomass_pred_cpue_strata'],
                      sd_log_pred = sdrep$sd[which(names(sdrep$value) == 'log_biomass_pred_cpue_strata')],
                      pred = exp(log_pred),
                      pred_lci = exp(log_pred - qnorm(1 - alpha_ci/2) * sd_log_pred),
                      pred_uci = exp(log_pred + qnorm(1 - alpha_ci/2) * sd_log_pred))
    }
  }

  # Model estimates of summed total biomass across strata and summed total cpue across strata when
  # available and appropriate to sum
  total_predicted_biomass <- NULL

  if(ncol(data$biomass_obs) == 1) {

    total_predicted_biomass <- ts_biomass_strata %>%
      dplyr::select(model_name, variable, year, pred, pred_lci, pred_uci) %>%
      dplyr::mutate(variable = 'tot_biomass_pred')

  } else {

    total_predicted_biomass <- tidyr::expand_grid(model_name = rema_model$input$model_name,
                                                  variable = 'tot_biomass_pred',
                                                  year = data$model_yrs) %>%
      dplyr::mutate(log_pred = sdrep$value[names(sdrep$value) == 'log_tot_biomass_pred'],
                    sd_log_pred = sdrep$sd[which(names(sdrep$value) == 'log_tot_biomass_pred')],
                    pred = exp(log_pred),
                    pred_lci = exp(log_pred - qnorm(1 - alpha_ci/2) * sd_log_pred),
                    pred_uci = exp(log_pred + qnorm(1 - alpha_ci/2) * sd_log_pred)) %>%
      dplyr::select(model_name, variable, year, pred, pred_lci, pred_uci)
  }

  # Prepare final output

  if(is.null(parameter_estimates)) {
    stop("Something went wrong... Please review the ?prepare_rema_input, ?fit_rema, and make sure REMA converged by running check_convergence(my_rema_model).")
  }

  if(!is.null(biomass_summary)) {
    biomass_by_strata <- biomass_summary
  } else {
    biomass_by_strata <- "Something went wrong... Please review the ?prepare_rema_input, ?fit_rema, and make sure REMA converged by running check_convergence(my_rema_model)."
  }

  if(!is.null(cpue_summary)) {
    cpue_by_strata <- cpue_summary
  } else {
    cpue_by_strata <- "REMA was fit only to biomass survey data, therefore no CPUE predictions available. If the user has a CPUE survey index and wants to fit to it, please see ?prepare_rema_input() for details."
  }

  if(!is.null(biomass_by_cpue_strata)) {
    biomass_by_cpue_strata <- biomass_by_cpue_strata
  } else {
    biomass_by_cpue_strata <- "'biomass_by_cpue_strata' is reserved for multi-survey scenarios when there are more biomass survey strata than CPUE survey strata, and the user wants predicted biomass at the same resolution as the CPUE survey index."
  }

  output <- list(parameter_estimates = parameter_estimates,
                 biomass_by_strata = biomass_by_strata,
                 cpue_by_strata = cpue_by_strata,
                 biomass_by_cpue_strata = biomass_by_cpue_strata,
                 total_predicted_biomass = total_predicted_biomass,
                 total_predicted_cpue = total_predicted_cpue)

  return(output)
}
