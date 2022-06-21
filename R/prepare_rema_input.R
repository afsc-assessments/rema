#' Prepare input data and parameters for REMA model
#'
#' After the data is read into R (either manually from a .csv or other data file
#' or by using \code{\link{read_admb_re}}, this function prepares the data and
#' parameter settings for \code{\link{fit_rema}}. The model can be set up to run
#' in single survey mode with one or more strata, or in multi-survey mode, which
#' uses an additional relative abundance index (i.e. cpue) to inform predicted
#' biomass. The optional inputs described below related to the CPUE survey data
#' or scaling parameter \code{q}, such as \code{cpue_dat} and \code{options_q}
#' are only used when \code{multi_survey = 1}. The function structure and
#' documentation is modeled after
#' \href{https://github.com/timjmiller/wham/blob/3a056359121bc1a911ed6a95c9203db4db456baa/R/prepare_wham_input.R}{wham::prepare_wham_input}.
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
#'     \code{pointer_PE_biomass = c(1, 1, 1)}. By default there is one unique
#'     log_PE estimated for each unique biomass survey stratum}
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
#'         number of PE parameters. For example, if three PE
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
#'     \item{$pointer_q_cpue}{An index to customize the assignment of q
#'     parameters to individual CPUE survey strata. Vector with length = number
#'     of CPUE strata, starting with an index of 1 and ending with the number of
#'     unique q parameters estimated. For example, if there are three CPUE
#'     survey strata and the user wanted to estimate only one q, they would
#'     specify \code{pointer_q_cpue = c(1, 1, 1)}. The recommended model
#'     configuration is to estimate one log_q for each CPUE survey stratum.}
#'     \item{$pointer_biomass_cpue_strata}{An index to customize the assignment
#'     of biomass predictions to individual CPUE survey strata. Vector with
#'     length = the number of biomass survey strata, starting with an index of 1
#'     and ending with the number of unique CPUE survey strata. This pointer
#'     only needs to be defined if the number of biomass and CPUE strata are not
#'     equal. The \code{pointer_biomass_cpue_strata} option allows the user to
#'     calculate predicted biomass at the CPUE survey strata level under the
#'     scenario where the biomass survey strata is at a higher resolution than
#'     the CPUE survey strata. For example, if there are 3 biomass survey strata
#'     that are represented by only 2 CPUE survey strata, the user may specify
#'     \code{pointer_biomass_cpue_strata = c(1, 1, 2)}. This specification would
#'     assign the first 2 biomass strata to the first CPUE strata, and the third
#'     biomass stratum to the second CPUE stratum. NOTE: there cannot be a
#'     scenario where there are more CPUE survey strata than biomass survey
#'     strata because the CPUE survey is used to inform the biomass survey
#'     trend. An error will be thrown if
#'     \code{q_options$pointer_biomass_cpue_strata} is not defined and the
#'     biomass and CPUE survey strata definitions are not the same.}
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
#' \code{zeros} allows the user to specify options for how to treat zero biomass
#' or CPUE survey observations. By default zero observations are treated as NAs
#' and a warning msg to that effect is returned to the console. \code{zeros}
#' allows the user to specify non-default zero assumptions using the following
#' list of entries:
#' \describe{
#'     \item{$assumption}{character, name of assumption using. Only three
#'     alternatives are currently implemented, \code{zeros = list(assumption =
#'     c("NA", "small_constant", "tweedie")}. \code{"NA"} is the default; this
#'     option assumes the zero estimates are failed surveys and removes them.
#'     \code{"small_constant"} is an ad hoc method where a fixed value is added
#'     to the zero with an assumed CV. By default, the small constant = 0.0001
#'     and the CV is the value entered by the user in the data. The user can
#'     change the assumed value and CV using \code{options_small_constant}.
#'     \code{"tweedie"} uses the Tweedie as the assumed error distribution of
#'     the survey data, which allows zeros. This alternative estimates
#'     additional power and dispersion parameters. The user can change initial
#'     values for these parameters or fix them at initial values using
#'     \code{options_tweedie}.}
#'     \item{$options_small_constant}{a vector length of
#'     two numeric values. The first value is the small constant to add to the
#'     zero observation, the second is the user-defined coefficient for this
#'     value. The user can specify the small value but use the input CV by
#'     specifying an NA for the second value. E.g., 'options_small_constant =
#'     c(0.0001, NA)'.}
#'     \item{$options_tweedie}{a list of entries to control initial or fixed
#'     values for Tweedie parameters. Currently, this argument accepts the
#'     following entries:}
#'         \describe{
#'         \item{$initial_pars}{In single-survey mode,
#'         \code{zeros$options_tweedie$initial_pars} must be a vector of numeric
#'         values with length = 2 \code{c(log_tweedie_dispersion,
#'         logit_tweedie_p)}. In multi-survey mode,
#'         \code{zeros$options_tweedie$initial_pars} must be a vector of numeric
#'         values with length = 4 \code{c(biomass survey log_tweedie_dispersion,
#'         logit_tweedie_p, cpue survey log_tweedie_dispersion,
#'         logit_tweedie_p)}. Initial values for log_tweedie_dispersion should
#'         be in log space. Initial values for logit_tweedie_p < -10 approach
#'         tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach
#'         tweedie_p = 2 (gamma).}
#'         \item{$fix_pars}{\code{zeros$options_tweedie$fix_pars} must be a
#'         vector of integer value(s) with the index value (starting at 1) of
#'         the log_tweedie_dispersion (1) or logit_tweedie_p (2) to be fixed.
#'         For example, in multi survey mode, if the user wants to fix both
#'         logit_tweedie_p's, they should enter \code{zeros = list(assumption =
#'         'tweedie', options_tweedie = list(fix_pars = c(2, 4)))}}
#'         }
#'     }
#'
#' @param admb_re list object returned from \code{\link{read_admb_re.R}}, which
#'   includes biomass survey data (\code{admb_re$biomass_dat}), optional cpue
#'   survey data (\code{admb_re$cpue_dat}), years for model predictions
#'   (\code{admb_re$model_yrs}), and model predictions of log biomass by strata
#'   in the correct format for input into REMA
#'   (\code{admb_re$init_log_biomass_pred}). If supplied, the user does not need
#'   enter biomass_dat or cpue_dat.
#' @param biomass_dat data.frame of biomass survey data in long format with the
#'   following columns:
#'   \describe{
#'   \item{\code{strata}}{character; the survey name, survey region, management
#'   unit, or depth strata. Note that the user must include this column even if
#'   there is only one survey strata}
#'   \item{\code{year}}{integer; survey year. Note that the user only needs to
#'   include years for which there are observations (i.e. there is no need to
#'   supply \code{NULL} or \code{NA} values for missing survey years)}
#'   \item{\code{biomass}}{numeric; the biomass estimate/observation (e.g.
#'   bottom trawl survey biomass in mt). By default, if \code{biomass == 0},
#'   this value will be treated as an NA (i.e., a failed survey). If the user
#'   wants to make other assumptions about zeros (e.g. adding a small constant),
#'   they must define it in the data manually.}
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
#'   survey cpue or relative population number); By default, if \code{cpue == 0},
#'   this value will be treated as an NA (i.e., a failed survey). If the user
#'   wants to make other assumptions about zeros (e.g. adding a small constant),
#'   they must define it in the data manually.}
#'   \item{\code{cv}}{numeric; the coefficient of variation (CV) of the cpue
#'   estimate (i.e. sd(cpue)/cpue)}
#'   }
#' @param start_year (optional) integer value specifying the start year for
#'   estimation in the model; if \code{admb_re} is supplied, this value defaults
#'   to \code{start_year = min(admb_re$model_yrs)}; if \code{admb_re} is not
#'   supplied, this value defaults to the first year in either
#'   \code{biomass_dat} or \code{cpue_dat}
#' @param end_year (optional) integer value specifying the last year for
#'   estimation in the model; if \code{admb_re} is supplied, this value defaults
#'   to \code{end_year = max(admb_re$model_yrs)}; if \code{admb_re} is not
#'   supplied, this value defaults to the last year in either \code{biomass_dat}
#'   or \code{cpue_dat}
#' @param sum_cpue_index T/F, is the CPUE survey index able to be summed across
#'   strata to get a total CPUE survey index? For example, Longline survey
#'   relative population numbers (RPNs) are summable but longline survey numbers
#'   per hachi (CPUE) are not. Default = \code{FALSE}.
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
#' @param zeros (optional) define assumptions about how to treat zero biomass or
#'   CPUE observations (see details).
#'
#' @return This function returns a named list with the following components:
#'   \describe{
#'   \item{\code{data}}{Named list of data, passed to
#'   \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'   \item{\code{par}}{Named list of parameters, passed to
#'   \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'   \item{\code{map}}{Named list defining how to optionally collect and fix
#'   parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'   \item{\code{random}}{Character vector of parameters to treat as random
#'   effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'   \item{\code{model_name}}{Name of stock or other identifier for REMA model}
#'   \item{\code{biomass_dat}}{A tidied long format data.frame of the biomass
#'   survey observations and associated CVs by strata. This data.frame will be
#'   'complete' in that it will include all modeled years, with missing values
#'   treated as NAs. Note that this data.frame could differ from the
#'   \code{admb_re$biomass_dat} or input \code{biomass} if assumptions about
#'   zero biomass observations are different between the ADMB model and what the
#'   user specifies for REMA. The user can change their assumptions about zeros
#'   using the \code{zeros} argument.}
#'   \item{\code{cpue_dat}}{If optional CPUE survey data are provided and
#'   \code{multi_survey = 1}, this will be a tidied long-format data.frame of
#'   the CPUE survey observations and associated CVs by strata. This data.frame
#'   will be 'complete' in that it will include all modeled years, with missing
#'   values treated as NAs. Note that this data.frame could differ from the
#'   \code{admb_re$biomass_dat} or input \code{biomass} if assumptions about
#'   zero CPUE observations are different between the ADMB model and what the
#'   user specifies for REMA. The user can change their assumptions about zeros
#'   using the \code{zeros} argument. If optional CPUE survey data are not
#'   provided or \code{multi_survey = 0}, this object will be \code{NULL}.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # place holder for example code
#' }
prepare_rema_input <- function(model_name = 'REMA for unnamed stock',
                               multi_survey = 0,
                               admb_re = NULL,
                               biomass_dat = NULL,
                               cpue_dat = NULL,
                               sum_cpue_index = FALSE,
                               start_year = NULL,
                               end_year = NULL,
                               wt_biomass = NULL,
                               wt_cpue = NULL,
                               PE_options = NULL,
                               q_options = NULL,
                               zeros = NULL) {

  # model_name = 'REMA for unnamed stock'
  # multi_survey = 0 # should this be 0 or FALSE instead of null?
  # admb_re = admb_re
  # # biomass_dat = NULL
  # # cpue_dat = NULL
  # sum_cpue_index = FALSE
  # start_year = NULL
  # end_year = NULL
  # wt_biomass = NULL
  # wt_cpue = NULL
  # PE_options = NULL
  # q_options = NULL
  # zeros = NULL

  data = list()
  par = list()
  map = list()
  random = character()

  input = list(data = data,
               par = par,
               map = map,
               random = random,
               model_name = model_name)

  # Fitting to more than one survey? default = no
  if(is.null(multi_survey)) {
    input$data$multi_survey = 0
  } else if (isTRUE(multi_survey) | multi_survey == 1) {
    input$data$multi_survey = 1
  } else {
    input$data$multi_survey = 0
  }

  # biomass and cpue survey data
  if(!is.null(admb_re)) {
    biomass_dat <- admb_re$biomass_dat
    cpue_dat <- admb_re$cpue_dat
  }

  # model years (years for predictions)
  if(!is.null(admb_re)) {

    model_yrs <- admb_re$model_yrs
    input$data$model_yrs <- model_yrs

  } else if(!is.null(cpue_dat)) {

    if(is.null(start_year)) start_year <- min(c(min(biomass_dat$year), min(cpue_dat$year)))
    if(is.null(end_year)) end_year <- max(c(max(biomass_dat$year), max(cpue_dat$year)))
    model_yrs <- start_year:end_year
    input$data$model_yrs <- model_yrs

  } else {
    if(is.null(start_year)) start_year <- min(biomass_dat$year)
    if(is.null(end_year)) end_year <- max(biomass_dat$year)
    model_yrs <- start_year:end_year
    input$data$model_yrs <- model_yrs
  }

  # see set_zeros.R
  # zeros <- set_zeros_default(zeros, dat = biomass_dat)
  if(is.null(zeros) & any(biomass_dat$biomass == 0, na.rm = TRUE)) {
    warning("The user has entered a zero observation for the survey data but has not explicitly defined an assumption for treatment of zeros in the model. By default, this observation will be removed (i.e. treated as an NA or failed survey). If the user wants to make another assumption (e.g. add a small constant or use the Tweedie distribution to model observation error), they can do so by defining an assumption in the 'zeros' argument. See details in ?prepare_rema_input() for more information.")
  }

  if(is.null(zeros$assumption)) {
    zeros <- list()
    zeros$assumption = 'NA'
  }

  biomass_dat = set_zeros(zeros = zeros,
                          dat = biomass_dat)

  # remove NAs, expand biomass survey data, and create biomass input for TMB
  biomass_dat <- biomass_dat %>%
    dplyr::filter(!is.na(biomass))

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

  # CPUE survey observations
  if((input$data$multi_survey == 1) & !is.null(cpue_dat)) {

    # zero survey cpue specs
    cpue_dat <- set_zeros(zeros = zeros,
                          dat = cpue_dat)

    # remove NAs
    cpue_dat <- cpue_dat %>%
      dplyr::filter(!is.na(cpue))

    cpue <- cpue_dat %>%
      tidyr::expand(year = model_yrs, strata) %>%
      dplyr::left_join(cpue_dat)

  } else if ((input$data$multi_survey == 1) & is.null(cpue_dat)){
    stop(paste("user defined multi_survey as TRUE but did not provide CPUE survey data in the admb_re list or as a dataframe in the cpue_dat argument."))

  } else {

    cpue <- expand.grid(year = model_yrs,
                        strata = unique(biomass_dat$strata),
                        cpue = NA,
                        cv = NA)
  }

  if((input$data$multi_survey == 0) & !is.null(cpue_dat)) {
    warning(paste('user defined multi_survey as FALSE but provided CPUE survey data. REMA will run in single survey (i.e. RE) mode and will not fit to the CPUE data. change to multi_survey = TRUE if you want to fit to survey CPUE data.'))
  }

  cpue_input <- cpue %>%
    tidyr::pivot_wider(id_cols = c("year"), names_from = "strata",
                       values_from = "cpue", values_fill = NA) %>%
    dplyr::mutate(value = "mu") %>%
    dplyr::bind_rows(cpue %>%
                       tidyr::pivot_wider(id_cols = c("year"), names_from = "strata",
                                          values_from = "cv", values_fill = NA) %>%
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

  # cpue index summable across strata?
  if(isTRUE(sum_cpue_index)){
    input$data$sum_cpue_index <- 1
  } else {
    input$data$sum_cpue_index <- 0 # default
  }

  # define default values for remaining data, par, and map list objects for TMB
  input <- set_defaults(input, admb_re = admb_re)

  # user-defined process error (PE) options
  input <- set_PE_options(input, PE_options)

  # user-defined scaling parameter (q) options
  input <- set_q_options(input, q_options)

  # set tweedie starting values, user-defined options
  input <- set_tweedie(input, zeros)

  if(!is.null(wt_biomass)) {
    input$data$wt_biomass <- wt_biomass
  }

  if(!is.null(wt_cpue)) {
    input$data$wt_cpue <- wt_cpue
  }

  # output tidied version of the biomass and cpue data
  input$biomass_dat <- biom %>%
    dplyr::arrange(strata, year)
  input$cpue_dat <- cpue %>%
    dplyr::arrange(strata, year)
  if(input$data$multi_survey == 0) {
    input$cpue_dat <- NULL
  }

  if(length(input$data$model_yrs) != nrow(input$par$log_biomass_pred)) {
    stop(paste0("Incorrect model dimensions! The number of rows for the log_biomass_pred input (nrow = ", nrow(input$par$log_biomass_pred), ") starting values (i.e. 'biomsd' in the rwout.rep file) and length of model years (", length(input$data$model_yrs), ") must match. The user can adjust the model 'start_year' and 'end_year' as needed."))
  }

  tstzeros <- as.vector(input$data$biomass_cv)
  if(any(tstzeros[!is.na(tstzeros)] <= 0)) {
    stop('The user has supplied a CV <= 0 for at least one biomass survey estimate, which is an invalid input.')
  }
  tstzeros <- as.vector(input$data$cpue_cv)
  if(any(tstzeros[!is.na(tstzeros)] <= 0)) {
    stop('The user has supplied a CV <= 0 for at least one CPUE survey observation, which is an invalid input.')
  }
  tstnegs <- as.vector(input$data$biomass_obs)
  if(any(tstnegs[!is.na(tstnegs)] < 0)) {
    stop('The user has supplied a negative value for at least one biomass survey estimate, which is an invalid input.')
  }
  tstnegs <- as.vector(input$data$cpue_obs)
  if(any(tstnegs[!is.na(tstnegs)] < 0)) {
    stop('The user has supplied a negative value for at least one biomass survey estimate, which is an invalid input.')
  }
  return(input)
}

