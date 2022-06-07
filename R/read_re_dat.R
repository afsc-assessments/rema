#' Convert ADMB version of the RE model output to REMA inputs
#'
#' Read the report file from the ADMB version of the RE model (rwout.rep) and
#' convert it into long format survey data estimates with CVs for input into
#' REMA.
#'
#' @param filename name of ADMB output file to be read (e.g. rwout.rep)
#' @param biomass_strata_names (optional) a vector of character names
#'   corresponding to the names of the biomass survey strata. Vector should be
#'   in the same order as the columns of srv_est in rwout.rep
#' @param cpue_strata_names (optional) a vector of character names corresponding
#'   to the names of the CPUE survey strata. Vector should be in the same order
#'   as the columns of srv_est_LL in rwout.rep in the version of the ADMB RE
#'   model that accepts an additional survey index
#' @return object of type "list" with biomass optional cpue survey data in
#'   long format, and initial parameter values for log_biomass_pred (the random
#'   effects matrix), ready for input into REMA
#' @export
read_re_dat <- function(filename,
                        biomass_strata_names = NULL,
                        cpue_strata_names = NULL) {

  # ex <- 'goasst.rep' # rema example
  # ex <- 'bsaisst.rep' # rem example
  ex <- 'aisr.rep' # re example
  # biomass_strata_names <- c('CGOA1', 'CGOA2', 'CGOA3',
  #                           'EGOA1', 'EGOA2', 'EGOA3',
  #                           'WGOA1', 'WGOA2', 'WGOA3')
  # cpue_strata_names <- c('CGOA', 'EGOA', 'WGOA')
  # filename <- re_ex

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

  # for multivariate models, assign strata names and check dimensions
  if(re_version %in% c('rem', 'rema')) {

    if(is.null(biomass_strata_names)) { biomass_strata_names <- paste('biomass_strata', 1:ncol(x$srv_est), sep = '_')
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
    biomass_est <- data.frame(strata = 'biomass_strata_1', year = x$yrs_srv, biomass = x$srv_est)
    biomass_cv <- data.frame(strata = 'biomass_strata_1', year = x$yrs_srv, cv = x$srv_sd)
  }

  # check for NAs, unique srv_yrs, number of strata
  biomass_dat <- biomass_est %>%
    dplyr::left_join(biomass_cv) %>%
    dplyr::arrange(strata, year) %>%
    dplyr::mutate(cv = ifelse(cv > 5, cv / biomass, cv)) %>%
    dplyr::select(strata, year, biomass, cv)

  if(re_version == 'rema') {

    if(is.null(cpue_strata_names)) { cpue_strata_names <- paste('cpue_strata', 1:ncol(x$srv_est_LL), sep = '_')
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
      dplyr::mutate(cv = ifelse(cv > 5, cv / cpue, cv)) %>%
      dplyr::select(strata, year, cpue, cv)

  } else {
    cpue_dat <- NULL
  }

  # initial values for log_biomass_pred
  init_log_biomass_pred <- x$biomsd
  if(is.vector(init_log_biomass_pred)) {
    init_log_biomass_pred <- as.matrix(init_log_biomass_pred)
  }
  re_dat <- list(biomass_dat = biomass_dat,
                 cpue_dat = cpue_dat,
                 init_log_biomass_pred = init_log_biomass_pred)

  return(re_dat)
}
