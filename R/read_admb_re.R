#' Convert ADMB version of the RE model data and output to REMA inputs
#'
#' Read the report file from the ADMB version of the RE model (rwout.rep) and
#' convert it into long format survey data estimates with CVs for input into
#' REMA.
#'
#' @param filename name of ADMB output file to be read (e.g. rwout.rep)
#' @param model_name (optional) Name of stock and identifier for the ADMB
#'   version of the RE model. Defaults to 'ADMB RE'
#' @param biomass_strata_names (optional) a vector of character names
#'   corresponding to the names of the biomass survey strata. Vector should be
#'   in the same order as the columns of srv_est in rwout.rep
#' @param cpue_strata_names (optional) a vector of character names corresponding
#'   to the names of the CPUE survey strata. Vector should be in the same order
#'   as the columns of srv_est_LL in rwout.rep in the version of the ADMB RE
#'   model that accepts an additional survey index
#'
#' @return object of type "list" with biomass optional cpue survey data in long
#'   format, and initial parameter values for log_biomass_pred (the random
#'   effects matrix), ready for input into REMA

#' @return a list with the following items:
#'   \describe{
#'     \item{\code{$biomass_dat}}{A dataframe of biomass survey data with
#'     strata, year, biomass estimates, and CVs. Note that the CVs have been
#'     back-transformed to natural space.}
#'     \item{\code{$cpue_dat}}{Optional dataframe of CPUE survey data with
#'     strata, year, CPUE estimates, and CVs. Note that the CVs have been
#'     back-transformed to natural space.}
#'     \item{\code{$model_yrs}}{Vector of prediction years.}
#'     \item{\code{$init_log_biomass_pred}}{Matrix of initial parameter values for
#'     log_biomass_pred (the random effects matrix), ready for input into REMA.}
#'     \item{\code{$admb_re_results}}{A list of ADMB RE model results ready for
#'     comparison with REMA models using compare_rema_models(). User beware...
#'     there are many, many versions of the RE.tpl in existence and individual
#'     variances may cause errors in this output.}
#'     }
#' @export
#'
#' @examples
#' \dontrun{
#' # place holder for example code
#' }
read_admb_re <- function(filename,
                         model_name = 'Unnamed ADMB RE model',
                         biomass_strata_names = NULL,
                         cpue_strata_names = NULL) {

  # filename <- 'inst/example_data/goasst_rwout.rep' # rema example 9 biomass strata, 3 cpue strata
  # filename <- 'inst/example_data/goasr_rwout.rep' # rema example 3 biomass strata, 3 cpue strata
  # filename <- 'inst/example_data/bsaisst_rwout.rep' # rem example
  # filename <- 'inst/example_data/aisr_rwout.rep' # re example

  x <- read_rep(fn = filename)

  # assign model version
  if(any(grepl('LL',names(x)))) {
    re_version <- 'rema'
  } else if(is.matrix(x$srv_est)) {
    re_version <- 'rem'
  } else {
    re_version <- 're'
  }
  # re_version

  # check that the variables needed for prepare_rema_input() exist in the
  # rwout.rep file provided by the user
  rem_names <- c('yrs_srv', 'srv_est', 'srv_sd', 'biomsd', 'yrs')
  rema_names <- c('yrs_srv_LL', 'srv_est_LL', 'srv_sd_LL')

  missing_names <- NULL
  if(!all(rem_names %in% names(x))) {
    missing_names <- c(missing_names, rem_names[!rem_names %in% names(x)])
  }
  if(re_version == 'rema' & !all(rema_names %in% names(x))) {
    missing_names <- c(missing_names, rema_names[!rema_names %in% names(x)])
  }
  if(!is.null(missing_names)) {
    stop(paste0("The following variable(s) are missing from the rwout.rep file provided by the user: ", toString(missing_names), ". These variables must be added to the rwout report file in order for read_admb_re() to function properly. This can be achieved by adding another write_R() statement in the re.tpl file to include the missing variables to the rwout.rep file, or it could be added manually."))
  }

  # assign strata names and check dimensions
  if(re_version %in% c('rem', 'rema')) {

    if(is.null(biomass_strata_names)) {
      biomass_strata_names <- paste('biomass_strata', 1:ncol(x$srv_est), sep = '_')
    }

    if(length(biomass_strata_names) != ncol(x$srv_est)) {
      stop(paste("the number of 'biomass_strata_names' provided by the user does not match the number of columns for 'srv_est' in the rwout.rep. the dimensions must match."))
    }

    biomass_est <- cbind(data.frame(x$yrs_srv), x$srv_est)
    names(biomass_est) <- c('year', biomass_strata_names)
    biomass_est <- biomass_est %>%
      tidyr::pivot_longer(cols = -year, names_to = 'strata', values_to = 'biomass') %>%
      dplyr::mutate(biomass = ifelse(biomass == -9, NA, biomass))

    biomass_cv <- cbind(data.frame(x$yrs_srv), x$srv_sd)
    names(biomass_cv) <- c('year', biomass_strata_names)
    biomass_cv <- biomass_cv %>%
      tidyr::pivot_longer(cols = -year, names_to = 'strata', values_to = 'cv') %>%
      dplyr::mutate(cv = ifelse(cv == -9, NA, cv))

  } else {

    if(is.null(biomass_strata_names)) {
      biomass_est <- data.frame(strata = 'biomass_strata_1', year = x$yrs_srv, biomass = x$srv_est)
      biomass_cv <- data.frame(strata = 'biomass_strata_1', year = x$yrs_srv, cv = x$srv_sd)
    }

    if(length(biomass_strata_names) > 1) {
      stop(paste("the number of 'biomass_strata_names' provided by the user is greater than one but there is only one biomass survey (i.e. 'srv_est') in the rwout.rep. Please provide only one biomass survey stratum name."))
    }

    if(!is.null(biomass_strata_names)){
      biomass_est <- data.frame(strata = biomass_strata_names, year = x$yrs_srv, biomass = x$srv_est)
      biomass_cv <- data.frame(strata = biomass_strata_names, year = x$yrs_srv, cv = x$srv_sd)
    }
  }

  # check for NAs, unique srv_yrs, number of strata
  biomass_dat <- biomass_est %>%
    dplyr::left_join(biomass_cv) %>%
    dplyr::arrange(strata, year) %>%
    # transform cv back to normal space (they were transformed to
    # log_biomasss_sd inside the re.tpl)
    dplyr::mutate(cv = sqrt(exp(cv ^ 2) - 1)) %>%
    dplyr::select(strata, year, biomass, cv)

  if(re_version == 'rema') {

    if(is.null(cpue_strata_names)) {
      cpue_strata_names <- paste('cpue_strata', 1:ncol(x$srv_est_LL), sep = '_')
    }

    if(length(cpue_strata_names) != ncol(x$srv_est_LL)) {
      stop(paste("the number of 'cpue_strata_names' provided by the user does not match the number of columns for 'srv_est_LL' in the rwout.rep."))
    }

    cpue_est <- cbind(data.frame(x$yrs_srv_LL), x$srv_est_LL)
    names(cpue_est) <- c('year', cpue_strata_names)
    cpue_est <- cpue_est %>%
      tidyr::pivot_longer(cols = -year, names_to = 'strata', values_to = 'cpue') %>%
      dplyr::mutate(cpue = ifelse(cpue == -9, NA, cpue))

    cpue_cv <- cbind(data.frame(x$yrs_srv_LL), x$srv_sd_LL)
    names(cpue_cv) <- c('year', cpue_strata_names)
    cpue_cv <- cpue_cv %>%
      tidyr::pivot_longer(cols = -year, names_to = 'strata', values_to = 'cv') %>%
      dplyr::mutate(cv = ifelse(cv == -9, NA, cv))

    cpue_dat <- cpue_est %>%
      dplyr::left_join(cpue_cv) %>%
      dplyr::arrange(strata, year) %>%
      # transform cv back to normal space (they were transformed to
      # log_biomasss_sd inside the re.tpl)
      dplyr::mutate(cv = sqrt(exp(cv ^ 2) - 1)) %>%
      dplyr::select(strata, year, cpue, cv)

  } else {
    cpue_dat <- NULL
  }

  # initial values for log_biomass_pred
  init_log_biomass_pred <- x$biomsd
  if(is.vector(init_log_biomass_pred)) {
    init_log_biomass_pred <- as.matrix(init_log_biomass_pred)
  }

  # years for predictions
  model_yrs <- x$yrs

  # check that the variables needed for compare_rema_models() exist in the
  # rwout.rep file provided by the user
  re_names <- c('biomA', 'LCI', 'UCI') # biomass model predictions
  missing_names <- NULL
  if(!all(re_names %in% names(x))) {
    missing_names <- c(missing_names, re_names[!re_names %in% names(x)])
  }
  rem_names <- c('biom_TOT', 'biom_TOT_LCI', 'biom_TOT_UCI')
  if(re_version %in% c('rem', 'rema') & !all(rem_names %in% names(x))) {
    missing_names <- c(missing_names, rem_names[!rem_names %in% names(x)])
  }
  rema_names <- c('biom_TOT_LL', 'biom_TOT_LCI_LL', 'biom_TOT_UCI_LL') # CPUE model predictions
  if(re_version %in% c('rema') & !all(rema_names %in% names(x))) {
    missing_names <- c(missing_names, rema_names[!rema_names %in% names(x)])
  }
  if(!is.null(missing_names)) {
    warning(paste0("The following variable(s) are missing from the rwout.rep file provided by the user: ", toString(missing_names), ". Although these variables are not needed to run REMA using prepare_rema_input() and fit_rema(), they will be needed to compare the results of the ADMB version of the RE model to REMA using compare_rema_models(). The user can add these variables to the rwout.rep file by coding in additional write_R() statements to the existing re.tpl file."))
  }

  # admb_re_results

  biomass_by_strata <- NULL

  if(!is.null(x$biomA)){
    biomass_by_strata <- cbind(data.frame(model_name = model_name,
                                               variable = 'biomass_pred',
                                               year = x$yrs),
                                    x$biomA)
    names(biomass_by_strata) <- c('model_name', 'variable', 'year', unique(biomass_dat$strata))
    biomass_by_strata <- biomass_by_strata %>%
      tidyr::pivot_longer(cols = -c(model_name, variable, year), names_to = 'strata', values_to = 'pred') %>%
      dplyr::arrange(strata, year)

    # upper and lower model 95% CIs
    if(!is.null(x$LCI) & !is.null(x$UCI)) {

      tmplci <- cbind(data.frame(year = x$yrs), x$LCI)
      names(tmplci) <- c('year', unique(biomass_dat$strata))
      tmplci <- tmplci %>%
        tidyr::pivot_longer(cols = -c(year), names_to = 'strata', values_to = 'pred_lci') %>%
        dplyr::arrange(strata, year)

      tmpuci <- cbind(data.frame(year = x$yrs), x$UCI)
      names(tmpuci) <- c('year', unique(biomass_dat$strata))
      tmpuci <- tmpuci %>%
        tidyr::pivot_longer(cols = -c(year), names_to = 'strata', values_to = 'pred_uci') %>%
        dplyr::arrange(strata, year)

      biomass_by_strata <- biomass_by_strata %>%
        dplyr::left_join(tmplci) %>%
        dplyr::left_join(tmpuci) %>%
        dplyr::arrange(strata, year)
    }
  } else {
    biomass_by_strata <- "The rwout.rep file provided by the user did not have 'biomA', 'LCI', 'UCI', the analagous variables in the ADMB version of the RE model to biomass_by_strata."
  }

  total_predicted_biomass <- NULL

  if(re_version == 're' & !is.data.frame(biomass_by_strata)) {
    total_predicted_biomass <- biomass_by_strata

  } else if(re_version == 're' & is.data.frame(biomass_by_strata)) {
    total_predicted_biomass <- biomass_by_strata  %>%
      dplyr::select(-strata) %>%
      dplyr::mutate(variable = 'tot_biomass_pred',
                    pred_sd = NA) %>%
      dplyr::select(model_name, variable, year, pred, pred_sd, pred_lci, pred_uci)

  } else if(!is.null(x$biom_TOT) & !is.null(x$biom_TOT_LCI) & !is.null(x$biom_TOT_UCI)) {

    total_predicted_biomass <- data.frame(model_name = model_name,
                                          variable = 'tot_biomass_pred',
                                          year = x$yrs,
                                          pred = x$biom_TOT,
                                          pred_sd = NA,
                                          pred_lci = x$biom_TOT_LCI,
                                          pred_uci = x$biom_TOT_UCI)
  } else {
    total_predicted_biomass <- "The rwout.rep file provided by the user did not have 'biom_TOT', 'biom_TOT_LCI', or 'biom_TOT_UCI', the analagous variables in the ADMB version of the RE model to total_predicted_biomass.  Please check the rwout.rep file."
  }

  biomass_by_cpue_strata <- "Unfortunately 'biomass_by_cpue_strata' is not output by default from the ADMB version of the RE model and is not readily available for comparison to REMA models. Please contact the author(s) of this package for more information."
  cpue_by_strata <- "Unfortunately 'cpue_by_strata' is not output by default from the ADMB version of the RE model and is not readily available for comparison to REMA models. Please contact the author(s) of this package for more information."

  total_predicted_cpue <- NULL

  if(re_version != 'rema') {
    total_predicted_cpue <- "'total_predicted_cpue' is not applicable to models that are not fit to CPUE data."
  } else if(re_version == 'rema' & !is.null(x$biom_TOT_LL) & !is.null(x$biom_TOT_LCI_LL) & !is.null(x$biom_TOT_UCI_LL)){

    total_predicted_cpue <- data.frame(model_name = model_name,
                                          variable = 'tot_cpue_pred',
                                          year = x$yrs,
                                          pred = x$biom_TOT_LL,
                                          pred_sd = NA,
                                          pred_lci = x$biom_TOT_LCI_LL,
                                          pred_uci = x$biom_TOT_UCI_LL)

  } else {
    total_predicted_cpue <- "The rwout.rep file provided by the user did not have 'biom_TOT_LL', 'biom_TOT_LCI_LL', or 'biom_TOT_UCI_LL', the analagous variables in the ADMB version of the RE model to total_predicted_cpue. Please check the rwout.rep file."
  }

  # Join to the ADMB RE data
  if(is.data.frame(biomass_by_strata)){

    alpha_ci <- 0.05
    biomass_by_strata <- biomass_by_strata %>%
      dplyr::mutate(pred_sd = NA) %>%
      dplyr::select(model_name, strata, variable, year, pred, pred_sd, pred_lci, pred_uci) %>%
      dplyr::left_join(biomass_dat %>%
                         dplyr::rename(obs = biomass,
                                       obs_cv = cv) %>%
                         # most common assumption in RE.tpl is to treat zeros as NAs
                         dplyr::mutate(obs = ifelse(obs == 0, NA, obs)) %>%
                         dplyr::filter(!is.na(obs)) %>%
                         dplyr::mutate(log_obs = log(obs),
                                       sd_log_obs = sqrt(log(obs_cv ^ 2 + 1)),
                                       obs_lci = exp(log_obs - qnorm(1 - alpha_ci/2) * sd_log_obs),
                                       obs_uci = exp(log_obs + qnorm(1 - alpha_ci/2) * sd_log_obs)))
  }
  parameter_estimates <- "The rwout.rep file does not contain parameter estimates, therefore they are not readily available for comparison with REMA models. Please contact the author(s) of this package for more information."

  admb_re_results <- list(parameter_estimates = parameter_estimates,
                          biomass_by_strata = biomass_by_strata,
                          cpue_by_strata = cpue_by_strata,
                          biomass_by_cpue_strata = biomass_by_cpue_strata,
                          total_predicted_biomass = total_predicted_biomass,
                          total_predicted_cpue = total_predicted_cpue)

  admb_re <- list(biomass_dat = biomass_dat,
                  cpue_dat = cpue_dat,
                  model_yrs = model_yrs,
                  init_log_biomass_pred = init_log_biomass_pred,
                  admb_re_results = admb_re_results)

  return(admb_re)
}
