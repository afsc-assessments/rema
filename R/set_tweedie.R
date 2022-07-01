
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
#              options_small_constant = c(0.0001, 5), # assumed mean, assumed CV
#              options_tweedie = list(initial_pars = c(1, 1.6, 1, 1.6),
#                                     fix_pars = c(2))) #, 4))

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
    if(is.null(zeros$options_tweedie)) {

      if(data$multi_survey == 0) {
        # map$log_tweedie_dispersion <- as.factor(c(1))
        map$logit_tweedie_p <- as.factor(c(1))
      }

      if(data$multi_survey == 1) {
        # map$log_tweedie_dispersion <- as.factor(c(1, 2))
        map$logit_tweedie_p <- as.factor(c(1, 2))
      }
    }

    if(!is.null(zeros$options_tweedie)) {

      # user-defined changes to initial values: Up to two values in the
      # following order: biomass survey logit_tweedie_p, cpue survey logit_tweedie_p

      if(!is.null(zeros$options_tweedie$initial_pars)) {

        if(data$multi_survey == 0) {
          if(length(zeros$options_tweedie$initial_pars) != 1) stop("In single-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 1 [c(logit_tweedie_p)]. Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
          if(!any(is.numeric(zeros$options_tweedie$initial_pars))) stop("In single-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 1 [c(logit_tweedie_p)]. Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
          # par$log_tweedie_dispersion <- zeros$options_tweedie$initial_pars[1]
          par$logit_tweedie_p <- zeros$options_tweedie$initial_pars[1]
        }

        if(data$multi_survey == 1) {
          if(length(zeros$options_tweedie$initial_pars) != 2) stop("In multi-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 2 [c(biomass survey logit_tweedie_p, cpue survey logit_tweedie_p). Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
          if(!any(is.numeric(zeros$options_tweedie$initial_pars))) stop("In multi-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 2 [c(biomass survey logit_tweedie_p, cpue survey logit_tweedie_p). Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
          # par$log_tweedie_dispersion <- zeros$options_tweedie$initial_pars[c(1, 3)]
          # par$logit_tweedie_p <- zeros$options_tweedie$initial_pars[c(2, 4)]
          # par$logit_tweedie_p <- zeros$options_tweedie$initial_pars
        }
      }

      # user-defined fixed parameters: Up to two values in the
      # following order: biomass survey
      # logit_tweedie_p, cpue survey logit_tweedie_p
      if(!is.null(zeros$options_tweedie$fix_pars)) {

        # fix_dispersion <- zeros$options_tweedie$fix_pars[which(zeros$options_tweedie$fix_pars %% 2 == 1)]
        # fix_p <- zeros$options_tweedie$fix_pars[which(zeros$options_tweedie$fix_pars %% 2 == 0)]

        # single survey mode
        if(data$multi_survey == 0) {

          # checks
          if(any(!zeros$options_tweedie$fix_pars %in% c(1))) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tweedie_p (1) to be fixed. See ?prepare_rema_input for details.")
          if(length(zeros$options_tweedie$fix_pars) > 1) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tweedie_p (1) to be fixed. See ?prepare_rema_input for details.")
          # zeros$options_tweedie$fix_pars <- as.integer(zeros$options_tweedie$fix_pars)

          # search through fix_par inputs and create map
          map$logit_tweedie_p <- as.factor(c(1))

        }

        # two survey mode
        if(data$multi_survey == 1) {

          # checks
          if(any(!zeros$options_tweedie$fix_pars %in% c(1:2))) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tweedie_p (1 or 2 for biomass and cpue survey parameters) to be fixed. See ?prepare_rema_input for details.")
          if(length(zeros$options_tweedie$fix_pars) > 2) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the logit_tweedie_p (1 or 2 for biomass and cpue survey parameters) to be fixed. See ?prepare_rema_input for details.")

          # search through fix_par and create map... would love to use a more
          # elegant way to do this.

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


# set_tweedie <- function(input, zeros) {
#
#   data = input$data
#   par = input$par
#   map = input$map
#
#   if(zeros$assumption != 'tweedie' & !is.null(zeros$options_tweedie)) {
#     stop("User has defined options_tweedie but did not define assumption = 'tweedie'. Please use assumption = 'tweedie' in order to use the Tweedie distribution for the observation error in the model, or remove the options_tweedie.")
#   }
#
#   if(zeros$assumption == 'tweedie') {
#
#     data$obs_error_type = 1
#
#     # default tweedie - note that parameter initial values for tweedie were
#     # already defined in set_defaults.R
#     if(is.null(zeros$options_tweedie)) {
#
#       if(data$multi_survey == 0) {
#         map$log_tweedie_dispersion <- as.factor(c(1))
#         map$logit_tweedie_p <- as.factor(c(1))
#       }
#
#       if(data$multi_survey == 1) {
#         map$log_tweedie_dispersion <- as.factor(c(1, 2))
#         map$logit_tweedie_p <- as.factor(c(1, 2))
#       }
#     }
#
#     if(!is.null(zeros$options_tweedie)) {
#
#       # user-defined changes to initial values: Up to four values in the
#       # following order: biomass survey log_tweedie_dispersion,
#       # logit_tweedie_p, cpue survey log_tweedie_dispersion, logit_tweedie_p
#
#       if(!is.null(zeros$options_tweedie$initial_pars)) {
#
#         if(data$multi_survey == 0) {
#           if(length(zeros$options_tweedie$initial_pars) != 2) stop("In single-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 2 [c(log_tweedie_dispersion, logit_tweedie_p ). Initial values for log_tweedie_dispersion should be in log space. Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
#           if(!any(is.numeric(zeros$options_tweedie$initial_pars))) stop("In single-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 2 [c(log_tweedie_dispersion, logit_tweedie_p ). Initial values for log_tweedie_dispersion should be in log space. Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
#           par$log_tweedie_dispersion <- zeros$options_tweedie$initial_pars[1]
#           par$logit_tweedie_p <- zeros$options_tweedie$initial_pars[2]
#         }
#
#         if(data$multi_survey == 1) {
#           if(length(zeros$options_tweedie$initial_pars) != 4) stop("In multi-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 4 [c(biomass survey log_tweedie_dispersion, logit_tweedie_p, cpue survey log_tweedie_dispersion, logit_tweedie_p). Initial values for log_tweedie_dispersion should be in log space. Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
#           if(!any(is.numeric(zeros$options_tweedie$initial_pars))) stop("In multi-survey mode, zeros$options_tweedie$initial_pars must be a vector of numeric values with length = 4 [c(biomass survey log_tweedie_dispersion, logit_tweedie_p, cpue survey log_tweedie_dispersion, logit_tweedie_p). Initial values for log_tweedie_dispersion should be in log space. Initial values for logit_tweedie_p < -10 approach tweedie_p = 1 (zero-inflated Poisson), logit_tweedie_p > 10 approach tweedie_p = 2 (gamma). See ?prepare_rema_input for details.")
#           par$log_tweedie_dispersion <- zeros$options_tweedie$initial_pars[c(1, 3)]
#           par$logit_tweedie_p <- zeros$options_tweedie$initial_pars[c(2, 4)]
#         }
#       }
#
#       # user-defined fixed parameters: Up to four values in the
#       # following order: biomass survey log_tweedie_dispersion,
#       # logit_tweedie_p, cpue survey log_tweedie_dispersion, logit_tweedie_p
#       if(!is.null(zeros$options_tweedie$fix_pars)) {
#
#         fix_dispersion <- zeros$options_tweedie$fix_pars[which(zeros$options_tweedie$fix_pars %% 2 == 1)]
#         fix_p <- zeros$options_tweedie$fix_pars[which(zeros$options_tweedie$fix_pars %% 2 == 0)]
#
#         # single survey mode
#         if(data$multi_survey == 0) {
#
#           # checks
#           if(any(!zeros$options_tweedie$fix_pars %in% c(1:2))) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the log_tweedie_dispersion (1) or logit_tweedie_p  (2) to be fixed. See ?prepare_rema_input for details.")
#           if(length(zeros$options_tweedie$fix_pars) > 2) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the log_tweedie_dispersion (1) or logit_tweedie_p  (2) to be fixed. See ?prepare_rema_input for details.")
#           # zeros$options_tweedie$fix_pars <- as.integer(zeros$options_tweedie$fix_pars)
#
#           # search through fix_par inputs and create map
#           if(length(fix_dispersion) == 0) {
#             map$log_tweedie_dispersion <- as.factor(c(1))
#           } else {
#             map$log_tweedie_dispersion <- as.factor(NA)
#           }
#
#           if(length(fix_p) == 0) {
#             map$logit_tweedie_p <- as.factor(c(1))
#           } else {
#             map$logit_tweedie_p <- as.factor(NA)
#           }
#         }
#
#         # two survey mode
#         if(data$multi_survey == 1) {
#
#           # checks
#           if(any(!zeros$options_tweedie$fix_pars %in% c(1:4))) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the log_tweedie_dispersion (1 or 3 for biomass and cpue survey parameters) or logit_tweedie_p (2 or 4 for biomass and cpue survey parameters) to be fixed. See ?prepare_rema_input for details.")
#           if(length(zeros$options_tweedie$fix_pars) > 4) stop("zeros$options_tweedie$fix_pars must be a vector of integer value(s) with the index value (starting at 1) of the log_tweedie_dispersion (1 or 3 for biomass and cpue survey parameters) or logit_tweedie_p (2 or 4 for biomass and cpue survey parameters) to be fixed. See ?prepare_rema_input for details.")
#           # zeros$options_tweedie$fix_pars <- as.integer(zeros$options_tweedie$fix_pars)
#
#           # search through fix_par and create map... would love to use a more
#           # elegant way to do this.
#
#           # dispersion parameters
#           if(length(fix_dispersion) == 0) {   # no dispersion parameters fixed
#             map$log_tweedie_dispersion <- as.factor(c(1, 2))
#
#           } else if(length(fix_dispersion) == 1) { # one dispersion parameters fixed
#
#             if(fix_dispersion == 1) { # biomass dispersion fixed
#               map$log_tweedie_dispersion <- as.factor(c(NA, 1))
#             }
#             if(fix_dispersion == 3) { # cpue dispersion fixed
#               map$log_tweedie_dispersion <- as.factor(c(1, NA))
#             }
#           } else { # both fixed
#             map$log_tweedie_dispersion <- as.factor(c(NA, NA))
#           }
#
#           # p parameters
#           if(length(fix_p) == 0) {   # no p parameters fixed
#             map$logit_tweedie_p <- as.factor(c(1, 2))
#
#           } else if(length(fix_p) == 1) { # one dispersion parameters fixed
#
#             if(fix_p == 2) { # biomass p fixed
#               map$logit_tweedie_p <- as.factor(c(NA, 1))
#             }
#             if(fix_p == 4) { # cpue p fixed
#               map$logit_tweedie_p <- as.factor(c(1, NA))
#             }
#           } else { # both fixed
#             map$logit_tweedie_p <- as.factor(c(NA, NA))
#           }
#         }
#       }
#     }
#   }
#
#   input$data = data
#   input$par = par
#   input$map = map
#
#   return(input)
#
# }
