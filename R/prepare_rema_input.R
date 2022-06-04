#' Prepare input data and parameters for REMA model
#'
#' After the data is read into R (either manually from a .csv or other data file
#' or by using \code{\link{read_re_dat.R}}, this function prepares the data and
#' parameter settings for \code{\link{fit_rema.R}}. The model can be set up to run
#' in single survey mode with one or more strata, or in multi-survey mode, which
#' uses an additional relative abundance index (i.e. cpue) to inform predicted
#' biomass. The inputs and options described below related to the cpue survey
#' data or scaling parameter $q$, such as \code{cpue_dat} and \code{options_q}
#' are only needed when \code{multi_survey = 1}.
#'
#' \code{PE_options} allows the user to specify options for process error (PE) parameters. If
#' \code{NULL}, default PE specifications are used: one PE parameter is
#' estimated for each biomass survey strata, initial values for \code{log_PE}
#' are set to 1, and no penalties or priors are added. The user can modify the
#' default \code{PE_options} using the following list of entries:
#' \describe{
#'     \item{$pointer_PE_biomass}{An index to customize the assignment of PE
#'     parameters to individual biomass strata. Vector with length = number of
#'     biomass strata, starting with an index of 1 and ending with the number of
#'     unique PE estimated. For example, if there are three biomass survey
#'     strata and the user wants to estimate only one PE, they would specify
#'     \code{pointer_PE_biomass = c(1, 1, 1). By default there is one unique
#'     log_PE estimated for each unique biomass survey stratum}}
#'     \item{$initial_pars}{A vector of initial values for log_PE. The
#'     default initial value for each log_PE is 1.}
#'     \item{$fix_pars}{Option to fix PE parameters, where the user specifies
#'     the index value of the PE parameter they would like to fix at the initial
#'     value. For example, if there are three biomass survey strata, and the
#'     user wants to fix the \code{log_PE} for the second stratum but estimate
#'     the \code{log_PE} for the first and third strata they would specify
#'     \code{fix_pars = c(2)} Note that this option is not recommended.}
#'     \item{$penalty_options}{Warning: the following options are experimental
#'     and not well-tested. Options for penalizing the PE likelihood or adding a
#'     prior on \code{log_PE} include the following:
#'     \describe{
#'         \item{"none"}{(default) no penalty or prior used}
#'         \item{"wt"}{a multiplier on the PE and random effects component of
#'         the negative log likelihood. For example, nll = wt * nll, where wt =
#'         1.5 is specified as a single value in the penalty_values argument}
#'        \item{"squared_penalty"}{As implemented in an earlier version of the
#'        RE.tpl, this penalty prevents the PE from shrinking to zero. For
#'        example, \code{nll = nll + (log_PE + squared_penalty)^2}, where
#'        \code{squared_penalty = 1.5}. A vector of \code{squared_penalty}
#'        values is specified for each PE in the \code{penalty_values} argument}
#'        \item{"normal_prior"}{Normal prior in log space, where \code{nll = nll -
#'        dnorm(log_PE, pmu_log_PE, psig_log_PE, 1)} and \code{pmu_log_PE} and
#'        \code{psig_log_PE} are specified for each PE parameter in the
#'        \code{penalty_values} argument}
#'        }
#'     }
#'     \item{penalty_values}{user-defined values for the \code{penalty_options}.
#'     Each penalty type will is entered as follows:
#'     \describe{
#'         \item{"none"}{(default) NULL For example, \code{penalty_values = NULL}}
#'         \item{"wt"}{a single numeric value. For example,
#'         \code{penalty_values = 1.5}}
#'         \item{"squared_penalty"}{a vector of numeric values with length =
#'         number of estimated PE parameters. For example, if three PE
#'         parameters are being estimated and the user wants them to have the
#'         same penalty for each one, they would use \code{penalty_values =
#'         c(1.5, 1.5, 1.5)}}
#'         \item{"normal_prior"}{a vector of paired values for each PE
#'         parameter, where each vector pair is the prior mean of log_PE
#'         \code{pmu_log_PE} and the associated standard deviation
#'         \code{psig_log_PE}. For example, if three PE parameters are being
#'         estimated and the user wants them to have the same normal prior of
#'         log_PE ~ N(1.0, 0.08), \code{penalty_values = c(c(1.0, 0.08), c(1.0,
#'         0.08), c(1.0, 0.08))}}
#'         }
#' }}
#'
#' \code{q_options} allows the user to specify options for the CPUE survey
#' scaling parameters (q). If \code{multi_survey = 0} (default), no q parameters
#' are estimated regardless of what the user defines in \code{q_options}.
#' \code{multi_survey = 0} and \code{q_options = NULL}, default q specifications
#' are used: one q parameter is estimated for each CPUE survey strata, biomass
#' and CPUE surveys are assumed to share strata definitions (i.e.,
#' \code{biomass_dat} and \code{cpue_dat} have the same number of columns and
#' the columns represent the same strata), initial values for \code{log_q} are
#' set to 1, and no penalties or priors are added. The user can modify the
#' default \code{q_options} using the following list of entries:
#' \describe{
#'     \item{$q_model}{Options for defining q parameters by CPUE survey strata, where:
#'     \describe{
#'         \item{"strata-specific"}{(default) one q for each CPUE strata}
#'         \item{"custom"}{user-defined q, where \code{pointer_q_cpue} must be
#'         specified}
#'      }
#'      }
#'     \item{$pointer_q_cpue}{An index to customize the assignment of q
#'     parameters to individual CPUE survey strata. Vector with length = number
#'     of CPUE strata, starting with an index of 1 and ending with the number of
#'     unique q parameters estimated. For example, if there are three CPUE
#'     survey strata and the user wanted to estimate only one q, they would
#'     specify \code{pointer_q_cpue = c(1, 1, 1)}}
#'     \item{$pointer_q_biomass}{An index to customize the assignment of q
#'     parameters to individual biomass survey strata. Vector with length =
#'     number of biomass survey strata, starting with an index of 1 and ending
#'     with the number of unique q parameters estimated. This pointer only needs
#'     to be defined if the number of biomass and CPUE strata are not equal. The
#'     \code{pointer_q_biomass} option allows the user to calculate predicted
#'     biomass at the CPUE survey strata level under the scenario where the
#'     biomass survey strata is at a higher resolution than the CPUE survey
#'     strata. For example, if there are 3 biomass survey strata that are
#'     represented by only 2 CPUE survey strata, the user may specify
#'     \code{pointer_q_biomass = c(1, 1, 2)}. This specification would assign
#'     the first 2 biomass strata to the first q, and the third biomass stratum
#'     to the second q. NOTE: there cannot be a scenario where there are more
#'     CPUE survey strata than biomass survey strata because the CPUE survey is
#'     used to inform the biomass survey trend.}
#'     \item{$initial_pars}{A vector of initial values for \code{log_q}. The
#'     default initial value for each log_q is 1.}
#'     \item{$fix_pars}{Option to fix q parameters, where
#'     the user specifies the index value of the q parameter they would like to
#'     fix at the initial value. For example, if there are three CPUE survey
#'     strata, and the user wants to fix the \code{log_q} for the second
#'     stratum but estimate the \code{log_q} for the first and third strata
#'     they would specify \code{fix_pars = c(2)}}
#'     \item{$penalty_options}{Options for penalizing the q likelihood or adding
#'     a prior on \code{log_q} include the following:
#'      \describe{
#'         \item{"none"}{(default) no penalty or prior used}
#'         \item{"normal_prior"}{Warning, experimental and not well-tested.
#'         Normal prior in log space, where \code{nll = nll - dnorm(log_q,
#'         pmu_log_q, psig_log_q, 1)} and \code{pmu_log_q} and \code{psig_log_q}
#'         are specified for each q parameter in the \code{penalty_values}
#'         argument}
#'        }
#'     }
#'     \item{penalty_values}{user-defined values for the \code{penalty_options}.
#'     Each penalty type will is entered as follows:
#'     \describe{
#'         \item{"none"}{(default) NULL For example, \code{penalty_values =
#'         NULL}}
#'         \item{"normal_prior"}{a vector of paired values for each q parameter,
#'         where each vector pair is the prior mean of log_q \code{pmu_log_q}
#'         and the associated standard deviation \code{psig_log_q}. For example,
#'         if 2 q parameters are being estimated and the user wants them to have
#'         the same normal prior of log_q ~ N(1.0, 0.05), \code{penalty_values =
#'         c(c(1.0, 0.05), c(1.0, 0.05))}}
#'         }
#' }}
#'
#' @param multi_survey logical; if equal to 1 (TRUE), the model will fit to an additional
#'   cpue survey index if provided in \code{cpue_dat}. Default = FALSE
#' @param re_dat list object returned from \code{\link{read_re_dat.R}}, which
#'   includes biomass_dat and optional cpue_dat in the correct format for input
#'   into REMA. If supplied, the user does not need enter biomass_dat or cpue_dat.
#' @param biomass_dat data.frame of biomass survey data in long format with the
#'   following columns:
#'   \describe{
#'   \item{\code{strata}}{character; the survey name, survey region, management
#'   unit, or depth strata. Note that the user must include this column even if
#'   there is only one survey strata}
#'   \item{\code{year}}{integer; survey year. Note that the user only needs to
#'   include years for which there are observations (i.e. there is no need to
#'   supply \code{NULL} or \code{NA} values for missing survey years)}
#'   \item{\code{biomass}}{numeric; the biomass estimate/observation (e.g. bottom trawl survey biomass in mt)}
#'   \item{\code{cv}}{numeric; the coefficient of variation (CV) of the biomass
#'   estimate (i.e. sd(biomass)/biomass)}
#'   }
#' @param cpue_dat (optional) data.frame of relative abundance index (i.e. cpue) data in
#'   long format with the following columns:
#'   \describe{
#'   \item{\code{strata}}{character; the survey name, survey region, management
#'   unit, or depth strata (note that the user must include this column even if
#'   there is only one survey strata)}
#'   \item{\code{year}}{integer; survey year. Note that the user only needs to
#'   include years for which there are observations (i.e. there is no need to
#'   supply \code{NULL} or \code{NA} values for missing survey years)}
#'   \item{\code{cpue}}{numeric; the cpue estimate/observation (e.g. longline
#'   survey cpue or relative population number)}
#'   \item{\code{cv}}{numeric; the coefficient of variation (CV) of the cpue
#'   estimate (i.e. sd(cpue)/cpue)}
#'   }
#' @param start_year (optional) integer value specifying the start year for
#'   estimation in the model; defaults to the first year of data in either
#'   \code{biomass_dat} or \code{cpue_dat}
#' @param end_year (optional) integer value specifying the last year for
#'   estimation in the model; defaults to the last year of data in either
#'   \code{biomass_dat} or \code{cpue_dat}
#' @param wt_biomass (optional) a multiplier on the biomass survey data
#'   component of the negative log likelihood. For example, \code{nll =
#'   wt_biomass * nll}. Defaults to \code{wt_biomass = 1}
#' @param wt_cpue (optional) a multiplier on the CPUE survey data
#'   component of the negative log likelihood. For example, \code{nll =
#'   wt_cpue * nll}. Defaults to \code{wt_cpue = 1}
#' @param PE_options (optional) customize implementation of process error (PE)
#'   parameters, including options to share PE across biomass survey strata,
#'   change starting values, fix parameters, and add penalties or priors (see
#'   details)
#' @param q_options (optional) customize implementation of scaling parameters
#'   (q), including options to define q by biomass or cpue survey cpue strata,
#'   change starting values, fix parameters, and add penalties or priors (see
#'   details). only used when \code{multi_survey = 1}
#'
#' @return This function returns a named list with the following components:
#'   \describe{ \item{\code{data}}{Named list of data, passed to
#'   \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}} \item{\code{par}}{Named list
#'   of parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'   \item{\code{map}}{Named list defining how to optionally collect and fix
#'   parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'   \item{\code{random}}{Character vector of parameters to treat as random
#'   effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'   \item{\code{years}}{Numeric vector of years to fit REMA model} }
#' @export
prepare_rema_input <- function(multi_survey = NULL,
                               re_dat = NULL,
                               biomass_dat = NULL,
                               cpue_dat = NULL,
                               start_year = NULL,
                               end_year = NULL,
                               wt_biomass = NULL,
                               wt_cpue = NULL,
                               PE_options = NULL,
                               q_options = NULL) {
  data = list()
  par = list()
  map = list()
  random = character()
  input = list(data = data,
               par = par,
               map = map,
               random = random)

  # Fitting to more than one survey? default = no
  if(is.null(multi_survey)) {
    input$data$multi_survey = 0
  }
  if(isTRUE(multi_survey)) {
    input$data$multi_survey = 1
  } else {
    input$data$multi_survey = 0
  }

  # biomass and cpue survey data
  if(!is.null(re_dat)) {
    biomass_dat <- re_dat$biomass_dat
    cpue_dat <- re_dat$cpue_dat
  }

  # model years
  if(is.null(start_year)) {
    start_year <- min(c(min(biomass_dat$year), min(cpue_dat$year)))
  }
  if(is.null(end_year)) {
    end_year <- max(c(max(biomass_dat$year), max(cpue_dat$year)))
  }
  model_yrs <- start_year:end_year
  input$data$model_yrs <- model_yrs

  # expand biomass and cpue survey data
  biom <- biomass_dat %>%
    tidyr::expand(year = model_yrs, strata) %>%
    dplyr::left_join(biomass_dat)

  biom_input <- biom %>%
    tidyr::pivot_wider(id_cols = c("year"), names_from = "strata",
                values_from = "biomass", values_fill = NA) %>%
    dplyr::mutate(value = "mu") %>%
    dplyr::bind_rows(biom %>%
                tidyr::pivot_wider(id_cols = c("year"), names_from = "strata",
                            values_from = "cv", values_fill = NA) %>%
                dplyr::mutate(value = "cv")) %>%
    dplyr::arrange(value, year)

  if(multi_survey == 1 & !is.null(cpue_dat)) {
    cpue <- cpue_dat %>%
      tidyr::expand(year = model_yrs, strata) %>%
      dplyr::left_join(cpue_dat)

  } else if (multi_survey == 1 & is.null(cpue_dat)){
    stop(paste("user defined multi_survey as TRUE but did not provide CPUE survey data in the re_dat list or as a dataframe in the cpue_dat argument."))

    } else {
    cpue <- expand.grid(year = model_yrs,
             strata = unique(biomass_dat$strata),
             cpue = NA,
             cv = NA)
    }

  if(multi_survey == 0 & !is.null(cpue_dat)) {
    warning(paste('user defined multi_survey as FALSE but provided CPUE survey data. REMA will run in single survey mode and will not fit to the CPUE data. change to multi_survey = TRUE if you want to fit to survey CPUE data.'))
  }

  cpue_input <- cpue %>%
    tidyr::pivot_wider(id_cols = c("year"), names_from = "strata",
                       values_from = "cpue", values_fill = NA) %>%
    dplyr::mutate(value = "mu") %>%
    dplyr::bind_rows(cpue %>%
                       tidyr::pivot_wider(id_cols = c("year"), names_from = "strata",
                                          values_from = "cv", values_fill = 0) %>%
                       dplyr::mutate(value = "cv")) %>%
    dplyr::arrange(value, year)

  # input biomass survey data
  input$data$biomass_obs <- biom_input %>%
    dplyr::filter(value == 'mu') %>%
    dplyr::select(-year, -value) %>%
    as.matrix()

  input$data$biomass_cv <- biom_input %>%
    dplyr::filter(value == 'cv') %>%
    dplyr::select(-year, -value) %>%
    as.matrix()

  # input cpue survey data
  input$data$cpue_obs <- cpue_input %>%
    dplyr::filter(value == 'mu') %>%
    dplyr::select(-year, -value) %>%
    as.matrix()

  input$data$cpue_cv <- cpue_input %>%
    dplyr::filter(value == 'cv') %>%
    dplyr::select(-year, -value) %>%
    as.matrix()

  # define default values for remaining data, par, and map list objects for TMB
  input <- set_defaults(input)

  # user-defined process error options, either to default values or settings specified by
  # user
  set_PE_options <- function(input, PE_options) {

    data = input$data
    par = input$par
    map = input$map

    # user defined index for PE estimation (e.g. there are 3 strata but user
    # only wants to estimate 1 PE, pointer_PE_biomass = c(1, 1, 1))
    if(!is.null(PE_options$pointer_PE_biomass)) {
      if(length(PE_options$pointer_PE_biomass) != ncol(data$biomass_obs)) stop("Length of PE_options$pointer_PE_biomass must equal the number of biomass survey strata (e.g. length(unique(re_dat$biomass_dat$strata))")
      PE_options$pointer_PE_biomass <- as.integer(PE_options$pointer_PE_biomass)
      if(!any(is.integer(PE_options$pointer_PE_biomass))) stop("PE_options$pointer_PE_biomass must be a vector of integer values starting at 1 with a vector length equal the number of biomass survey strata (e.g. length(unique(re_dat$biomass_dat$strata))")
      data$pointer_PE_biomass <- (PE_options$pointer_PE_biomass)-1 # TMB started indexing at 0
      par$log_PE <- rep(1, length(unique(data$pointer_PE_biomass)))
      map$log_PE <- fill_vals(par$log_PE, NA)
    }

    # user defined initial values for log_PE parameters
    if(!is.null(PE_options$initial_pars)) {
      if(length(PE_options$initial_pars) != length(par$log_PE)) stop("PE_options$initial_pars must be a vector with length equal to the number of PE parameters in the model. Input values must be numeric and greater than zero, because log_PE parameters are estimated in log space.")
      if(any(PE_options$initial_pars <= 0)) stop("PE_options$initial_pars must be a vector with length equal to the number of PE parameters in the model. Input values must be numeric and greater than zero, because log_PE parameters are estimated in log space.")
      if(!any(is.numeric(PE_options$initial_pars))) stop("PE_options$initial_pars must be a vector with length equal to the number of PE parameters in the model. Input values must be numeric and greater than zero, because log_PE parameters are estimated in log space.")
      par$log_PE <- PE_options$initial_pars
      map$log_PE <- fill_vals(par$log_PE, NA)
    }

    # user defined fixed values of log_PE. note that this is not recommended.
    if(!is.null(PE_options$fix_pars)) {
      warning("Are you sure you want to fix process error parameters? This is not recommended but may be useful for sensitivity analysis or troubleshooting. User can fix log_PE value to something other than the default of log_PE = 1 using the PE_options$initial_pars setting.")

      PE_options$fix_pars <- as.integer(PE_options$fix_pars)
      if(!any(is.integer(PE_options$fix_pars))) stop("PE_options$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the log_PE to be fixed.")

      # Next - figure out how to map off certain parameters in the log_PE vector.
      par$log_PE[PE_options$fix_pars]
    }
  }
  # test values - remove when fxn is complete
  PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                    initial_pars = c(1.1, 1.1, 1.1),
                    fix_pars = c(2, 3),
                    penalty_options = 'none',
                    penalty_values = NULL)

  pointer_q_biomass

  return(input)

}

