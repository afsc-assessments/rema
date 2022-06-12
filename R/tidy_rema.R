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
#'     \item{\code{$total_predicted_biomass}}{A tidy, long format data.frame of total
#'     model predicted biomass summed across all biomass survey strata. If only one stratum is
#'     used (i.e. the univariate RE), this will return a character string
#'     directing the user to \code{output$biomass_by_strata}.}
#'     \item{\code{$total_predicted_cpue}}{A tidy, long format data.frame of
#'     total model predicted CPUE summed across all CPUE survey strata. If only
#'     one stratum is used (i.e. the univariate RE), this will return a
#'     character string directing the user to \code{output$cpue_by_strata}. If
#'     The CPUE survey index provided was defined as not summable in
#'     prepare_rema_input(), an character string will be returned explaining how
#'     to change this using the 'sum_cpue_index' in \code{?prepare_rema_input}
#'     if appropriate.}
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

  # rema_model = m

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

    pe_pars <- data.frame(model_name = input$model_name,
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

  # scaling parameters when available
  q_pars <- NULL
  if(data$multi_survey == 1 & length(pars$log_q) > 0) {
    q_pars <- data.frame(model_name = input$model_name,
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

  parameter_estimates <- NULL

  if(!is.null(q_pars)) {
    parameter_estimates <- dplyr::bind_rows(pe_pars, q_pars)
  } else {
    parameter_estimates <- pe_pars
  }

  # Derived quantities:

  # Model estimates of summed total biomass across strata and summed total cpue across strata when
  # available and appropriate to sum
  ts_totals <- NULL
  ts_totals <- data.frame(model_name = input$model_name,
                          variable = c('tot_biomass_pred', 'tot_cpue_pred')) %>%
    dplyr::right_join(tidyr::expand_grid(variable = unique(names(sdrep$value)[grepl('tot_', names(sdrep$value))]),
                                         year = data$model_yrs) %>%
                        dplyr::mutate(pred = sdrep$value[grepl('tot_', names(sdrep$value))],
                                      pred_sd = sdrep$sd[which(grepl('tot_', names(sdrep$value)))]) %>%
                        dplyr::mutate(pred_lci = pred - qnorm(1 - alpha_ci/2) * pred_sd,
                                      pred_uci = pred + qnorm(1 - alpha_ci/2) * pred_sd))

  # Model estimates of biomass by strata
  ts_biomass_strata <- NULL
  ts_biomass_strata <- tidyr::expand_grid(model_name = input$model_name,
                                          strata = colnames(data$biomass_obs),
                                          variable = unique(names(sdrep$value)[!grepl(c('tot_|cpue'), names(sdrep$value))]),
                                          year = data$model_yrs) %>%
    dplyr::mutate(pred = sdrep$value[!grepl(c('tot_|cpue'), names(sdrep$value))],
                  pred_sd = sdrep$sd[which(!grepl(c('tot_|cpue'), names(sdrep$value)))]) %>%
    dplyr::mutate(pred_lci = pred - qnorm(1 - alpha_ci/2) * pred_sd,
                  pred_uci = pred + qnorm(1 - alpha_ci/2) * pred_sd) %>%
    dplyr::mutate(pred_lci = ifelse(pred_lci < 0, 0, pred_lci))

  biomass_summary <- NULL
  biomass_summary <- ts_biomass_strata %>%
    dplyr::filter(variable == 'biomass_pred') %>%
    dplyr::left_join(rema_model$input$biomass_dat %>%
                       dplyr::rename(obs = biomass, obs_cv = cv) %>%
                       dplyr::mutate(obs_sd = obs_cv * obs,
                                     obs_lci = obs - qnorm(1 - alpha_ci/2) * obs_sd,
                                     obs_uci = obs + qnorm(1 - alpha_ci/2) * obs_sd) %>%
                       dplyr::mutate(obs_lci = ifelse(obs_lci < 0, 0, obs_lci)))

  # Model estimates of cpue by strata when available
  ts_cpue_strata <- NULL
  cpue_summary <- NULL
  biomass_by_cpue_strata <- NULL

  if(data$multi_survey == 1){
    ts_cpue_strata <- tidyr::expand_grid(model_name = input$model_name,
                                         strata = colnames(data$cpue_obs),
                                         variable = unique(names(sdrep$value)[names(sdrep$value) %in% c('cpue_pred')]),
                                         year = data$model_yrs) %>%
      dplyr::mutate(pred = sdrep$value[names(sdrep$value) %in% c('cpue_pred')],
                    pred_sd = sdrep$sd[which(names(sdrep$value) %in% c('cpue_pred'))]) %>%
      dplyr::mutate(pred_lci = pred - qnorm(1 - alpha_ci/2) * pred_sd,
                    pred_uci = pred + qnorm(1 - alpha_ci/2) * pred_sd) %>%
      plyr::mutate(pred_lci = ifelse(pred_lci < 0, 0, pred_lci))

    cpue_summary <- ts_cpue_strata %>%
      dplyr::filter(variable == 'cpue_pred') %>%
      dplyr::left_join(rema_model$input$cpue_dat %>%
                         dplyr::rename(obs = cpue, obs_cv = cv) %>%
                         dplyr::mutate(obs_sd = obs_cv * obs,
                                       obs_lci = obs - qnorm(1 - alpha_ci/2) * obs_sd,
                                       obs_uci = obs + qnorm(1 - alpha_ci/2) * obs_sd) %>%
                         dplyr::mutate(obs_lci = ifelse(obs_lci < 0, 0, obs_lci)))

    # if there are more biomass strata than cpue strata, get
    # predicted biomass by cpue strata for comparison at the same strata level
    if(any(grepl('biomass_pred_cpue_strata', unique(names(sdrep$value))))) {

      biomass_by_cpue_strata <- tidyr::expand_grid(model_name = input$model_name,
                                           strata = colnames(data$cpue_obs),
                                           variable = unique(names(sdrep$value)[names(sdrep$value) %in% c('biomass_pred_cpue_strata')]),
                                           year = data$model_yrs) %>%
        dplyr::mutate(pred = sdrep$value[names(sdrep$value) %in% c('biomass_pred_cpue_strata')],
                      pred_sd = sdrep$sd[which(names(sdrep$value) %in% c('biomass_pred_cpue_strata'))]) %>%
        dplyr::mutate(pred_lci = pred - qnorm(1 - alpha_ci/2) * pred_sd,
                      pred_uci = pred + qnorm(1 - alpha_ci/2) * pred_sd) %>%
        plyr::mutate(pred_lci = ifelse(pred_lci < 0, 0, pred_lci))

      # join biomass and cpue strata definitions

      # strata_lkup <- data.frame(strata = unique(rema_model$input$biomass_dat$strata),
      #                           pointer = data$pointer_q_biomass + 1) %>%
      #   dplyr::left_join(data.frame(cpue_strata = unique(rema_model$input$cpue_dat$strata),
      #                               pointer = data$pointer_q_cpue + 1))

      biomass_by_cpue_strata <- biomass_by_cpue_strata %>%
        dplyr::filter(variable == 'biomass_pred_cpue_strata') # %>%
      # following code commented out because in most cases it is inappropriate
      # to sum the biomass observations because there could be a missing value
      # in one summed strata but not the other within a given year
      # dplyr::left_join(rema_model$input$biomass_dat %>%
      #                    dplyr::left_join(strata_lkup) %>%
      #                    dplyr::mutate(strata = cpue_strata,
      #                                  var = (cv * biomass)^2) %>%
      #                    dplyr::select(-cpue_strata) %>%
      #                    dplyr::group_by(strata, year) %>%
      #                    dplyr::summarise(obs = sum(biomass),
      #                                     var = sum(var)) %>%
      #                    dplyr::ungroup() %>%
      #                    dplyr::mutate(obs_cv = sqrt(var) / obs,
      #                                  obs_sd = sqrt(var),
      #                                  obs_lci = obs - qnorm(1 - alpha_ci/2) * obs_sd,
      #                                  obs_uci = obs + qnorm(1 - alpha_ci/2) * obs_sd) %>%
      #                    dplyr::mutate(obs_lci = ifelse(obs_lci < 0, 0, obs_lci)) %>%
      #                    dplyr::select(-var))
    }
  }

  # Prepare final output

  if(is.null(ts_totals)) {
    stop("Something went wrong... Please review the ?prepare_rema_input, ?fit_rema, and make sure REMA converged by running check_convergence(my_rema_model).")
  }

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

  total_predicted_biomass <- ts_totals %>%
    dplyr::filter(variable == 'tot_biomass_pred')

  # if(length(unique(biomass_by_strata$strata)) == 1) {
  #   message("Just an fyi: only one biomass survey stratum was fit in REMA, therefore the predicted values in output$biomass_by_strata and output$total_predicted_biomass will be the same.")
  # }

  if(data$multi_survey == 0) {
    total_predicted_cpue <- "REMA was fit only to biomass survey data, therefore no CPUE predictions available. If the user has a CPUE survey index and wants to fit to it, please see ?prepare_rema_input() for details."
  } else if (data$sum_cpue_index == 0) {
    total_predicted_cpue <- "The CPUE survey index provided was defined as not summable in prepare_rema_input(). If the CPUE index is summable (e.g. Relative Population Numbers), please select sum_cpue_index = TRUE in prepare_rema_input(). See ?prepare_rema_input() for more details."
  } else {
    total_predicted_cpue <- ts_totals %>%
      dplyr::filter(variable == 'tot_cpue_pred')

    # if(length(unique(cpue_by_strata$strata) == 1)) {
    #   message("Just an fyi: only one CPUE survey stratum was fit in REMA, therefore the predicted values in output$cpue_by_strata and output$total_predicted_cpue will be the same.")
    # }
  }

  output <- list(parameter_estimates = parameter_estimates,
                 biomass_by_strata = biomass_by_strata,
                 cpue_by_strata = cpue_by_strata,
                 biomass_by_cpue_strata = biomass_by_cpue_strata,
                 total_predicted_biomass = total_predicted_biomass,
                 total_predicted_cpue = total_predicted_cpue)

  return(output)
}
