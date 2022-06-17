
# function to set default zeros values (currently treat as NAs) replace zeros in the biomass or cpue _dat objects based on
# user-defined assumptions. see 'zeros' details in ?prepare_rema_inputs

set_zeros_default <- function(zeros) {

  # set default assumptions
  if(is.null(zeros) & any(dat$biomass == 0, na.rm = TRUE)) {
    warning("The user has entered a zero observation for the biomass survey data but has not explicitly defined an assumption for treatment of zeros in the model. By default, this observation will be removed (i.e. treated as an NA or failed survey). If the user wants to make another assumption (e.g. add a small constant or use the Tweedie distribution to model observation error), they can do so by defining an assumption in the 'zeros' argument. See details in ?prepare_rema_input() for more information.")
  }

  if(is.null(zeros$assumption)) {
    zeros <- list()
    zeros$assumption = 'NA'
  }

  return(zeros)
}

set_zeros <- function(dat  # either biomass_dat or cpue_dat
){

  # tests
  # dat = data.frame(year = 1990:1994, strata = 'strata', biomass = c(NA, 0, 0, 1, 1), cv = rep(0.1, 5))
  # zeros = list(assumption = 'tweedie', # c('NA', 'small_constant', 'tweedie'),
  #              options_small_constant = c(0.0001, NA), # replace mean, CV from data
  #              # options_small_constant = c(0.0001, 5), # replace mean, replace CV
  #              options_tweedie = list(initial_pars = c(1, 1.6, 1, 1.6), # biomass survey (phi, rho); cpue survey (phi, rho)
  #                                     fix_pars = NULL)
  # )
  # zeros = NULL

  # set optional inputs for adding a small constant. user can define the value (defaults to 1e-4)
  default_small_constant <- 1e-4

  if(zeros$assumption == 'small_constant') {

    if(is.null(zeros$options_small_constant)) {
      tmp_small_constant <- default_small_constant
    } else {
      if(length(zeros$options_small_constant) != 2 & !is.numeric(zeros$options_small_constant)) {
        stop("For 'options_small_constant' the user must enter a vector length of two numeric values. The first value is the small constant to add to the zero observation, the second is the user-defined coefficient for this value. The user can specify the small value but use the input CV by specifying an NA for the second value. E.g., 'options_small_constant = c(0.0001, NA)'. See ?prepare_rema_inputs for details.")
      } else {

        # user small constant for zeros
        if(is.na(zeros$options_small_constant[1])) {
          tmp_small_constant <- default_small_constant
        } else {
          tmp_small_constant <- zeros$options_small_constant[1]
        }

        # user cv for zeros
        if(is.na(zeros$options_small_constant[2])) {
          tmp_cv <- NULL
        } else {
          tmp_cv <- zeros$options_small_constant[2]
        }
      }
    }
  }

  # set biomass ----
  if(any(names(dat) %in% c('biomass'))) {

    if(c(zeros$assumption == 'NA')) {
      dat <- dat %>%
        dplyr::mutate(biomass = ifelse(biomass == 0, NA, biomass))
    }

    if(c(zeros$assumption == 'small_constant')) {
      dat <- dat %>%
        dplyr::mutate(cv = ifelse(biomass == 0 & !is.null(tmp_cv), tmp_cv, cv),
                      biomass = ifelse(biomass == 0, tmp_small_constant, biomass))
    }
  }

  # set cpue ----
  if(any(names(dat) %in% c('cpue'))) {

    if(c(zeros$assumption == 'NA')) {
      dat <- dat %>%
        dplyr::mutate(cpue = ifelse(cpue == 0, NA, cpue))
    }

    if(c(zeros$assumption == 'small_constant')) {
      dat <- dat %>%
        dplyr::mutate(cv = ifelse(cpue == 0 & !is.null(tmp_cv), tmp_cv, cv),
                      cpue = ifelse(cpue == 0, tmp_small_constant, cpue))
    }
  }
  # dat
  return(dat)
}
