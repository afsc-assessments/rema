
# function to define Tweedie distribution starting values and fixed parameters
# within prepare_rema_input()

# tweedie_p_init = 1.6 # based on Dave Warton
# tweedie_p_upper = 2 # gamma
# tweedie_p_lower = 1 # zero inf poisson
# logit_tweedie_p = log((tweedie_p_init - tweedie_p_lower) / (tweedie_p_upper - tweedie_p_init))
# logit_tweedie_p
# tweedie_p = Type(1.0) + ((Type(2.0) - Type(1.0) / (Type(1.0) + exp(-logit_tweedie_p))))
# tweedie_p = 1 + (1 / (1 + exp(-logit_tweedie_p)))
# tweedie_p

# zeros = list(assumption = 'tweedie', # c('NA', 'small_constant', 'tweedie'),
#              # options_small_constant = c(0.0001, 5), # assumed mean, assumed CV
#              options_tweedie = list(zeros_cv = 1.8,
#                                     initial_pars = c(1.6),
#                                     fix_pars = c(1))) #, 4))

set_tweedie <- function(input, zeros) {

  data = input$data
  par = input$par
  map = input$map

  if(zeros$assumption != 'tweedie' & !is.null(zeros$options_tweedie)) {
    stop("User has defined options_tweedie but did not define assumption = 'tweedie'. Please use assumption = 'tweedie' in order to use the Tweedie distribution for the observation error in the model, or remove the options_tweedie.")
  }

  if(zeros$assumption == 'tweedie') {

    data$obs_error_type = 1

    # default tweedie - note that parameter initial values for tweedie were
    # already defined in set_defaults.R
    if(data$multi_survey == 0) {
        map$logit_tweedie_p <- as.factor(c(1))
    }

    if(data$multi_survey == 1) {
        map$logit_tweedie_p <- as.factor(c(1, 2))
    }

    if(!is.null(zeros$options_tweedie)) {

      # user-defined changes to zeros_cv, the assumed CV of zero biomass or cpue
      # survey observations. default = 1.5
      if(!is.null(zeros$options_tweedie$zeros_cv)) {

        if(length(zeros$options_tweedie$zeros_cv) != 1) stop("zeros$options_tweedie$zeros_cv accepts a single positive, non-zero value as an input. This input changes the assumed value for zero biomass or cpue survey observations. See ?prepare_rema_input for details.")
        if(!any(is.numeric(zeros$options_tweedie$zeros_cv))) stop("zeros$options_tweedie$zeros_cv accepts a single positive, non-zero value as an input. This input changes the assumed value for zero biomass or cpue survey observations. See ?prepare_rema_input for details.")
        if(!any(zeros$options_tweedie$zeros_cv > 0)) stop("zeros$options_tweedie$zeros_cv accepts a single positive, non-zero value as an input. This input changes the assumed value for zero biomass or cpue survey observations. See ?prepare_rema_input for details.")

        data$biomass_cv[data$biomass_obs==0] <- zeros$options_tweedie$zeros_cv
        data$cpue_cv[data$cpue_obs==0] <- zeros$options_tweedie$zeros_cv
      }

      # user-defined changes to initial values: Up to two values in the
      # following order: biomass survey logit_tweedie_p, cpue survey logit_tweedie_p

      if(!is.null(zeros$options_tweedie$initial_pars)) {

        if(data$multi_survey == 0) {
          if(length(zeros$options_tweedie$initial_pars) != 1) stop("In single-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 1 [c(logit_tweedie_p)]. Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
          if(!any(is.numeric(zeros$options_tweedie$initial_pars))) stop("In single-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 1 [c(logit_tweedie_p)]. Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
          par$logit_tweedie_p <- zeros$options_tweedie$initial_pars
        }

        if(data$multi_survey == 1) {
          if(length(zeros$options_tweedie$initial_pars) != 2) stop("In multi-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 2 [c(biomass survey logit_tweedie_p, cpue survey logit_tweedie_p). Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
          if(!any(is.numeric(zeros$options_tweedie$initial_pars))) stop("In multi-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 2 [c(biomass survey logit_tweedie_p, cpue survey logit_tweedie_p). Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
          par$logit_tweedie_p <- zeros$options_tweedie$initial_pars
        }
      }

      # user-defined fixed parameters: Up to two values in the
      # following order: biomass survey
      # logit_tweedie_p, cpue survey logit_tweedie_p
      if(!is.null(zeros$options_tweedie$fix_pars)) {

        # single survey mode
        if(data$multi_survey == 0) {

          # checks
          if(any(!zeros$options_tweedie$fix_pars %in% c(1))) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tweedie_p (1) to be fixed. See ?prepare_rema_input for details.")
          if(length(zeros$options_tweedie$fix_pars) > 1) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tweedie_p (1) to be fixed. See ?prepare_rema_input for details.")
          # zeros$options_tweedie$fix_pars <- as.integer(zeros$options_tweedie$fix_pars)

          # search through fix_par inputs and create map
          map$logit_tweedie_p <- as.factor(NA)

        }

        # two survey mode
        if(data$multi_survey == 1) {

          # checks
          if(any(!zeros$options_tweedie$fix_pars %in% c(1:2))) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tweedie_p (1 or 2 for biomass and cpue survey parameters) to be fixed. See ?prepare_rema_input for details.")
          if(length(zeros$options_tweedie$fix_pars) > 2) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tweedie_p (1 or 2 for biomass and cpue survey parameters) to be fixed. See ?prepare_rema_input for details.")

          # p parameters
          if(length(zeros$options_tweedie$fix_pars) == 2) {
            map$logit_tweedie_p <- as.factor(c(NA, NA))
          }
          if(length(zeros$options_tweedie$fix_pars) == 1) {

            if(zeros$options_tweedie$fix_pars == 1) {
              map$logit_tweedie_p <- as.factor(c(NA, 1))
            }

            if(zeros$options_tweedie$fix_pars == 2) {
              map$logit_tweedie_p <- as.factor(c(1, NA))
            }

          }

        }
      }
    }
  }

  input$data = data
  input$par = par
  input$map = map

  return(input)

}
