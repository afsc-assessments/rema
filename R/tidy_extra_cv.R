#' Tidy estimates of extra biomass or CPUE index CV
#'
#' Takes list output from \code{\link{tidy_rema}} and returns the same list with
#' enhanced versions of the \code{biomass_by_strata} and \code{cpue_by_strata}
#' when appropriate. These enhanced dataframes include three new columns,
#' \code{tot_log_obs_cv}, \code{tot_obs_lci}, and \code{tot_obs_uci}, which
#' represent combined log-space standard error and associated confidence
#' intervals that include both assumed and estimated additional observation
#' error.
#'
#' @param tidy_rema list out output from \code{\link{tidy_rema}}, which includes
#'   model results but also inputs
#' @param save (optional) logical (T/F) save figures as \code{filetype} in
#'   \code{path}. Default = FALSE. NOT YET IMPLEMENTED.
#' @param path (optional) directory path to location where figure files are to
#'   be saved if \code{save = TRUE}. NOT YET IMPLEMENTED.
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
#'     predicted and observed biomass by biomass survey strata. This data.frame
#'     is now enhanced with new columns that include log-space standard error
#'     and associated confidence intervals that account for additional estimated
#'     observation error.}
#'     \item{\code{$cpue_by_strata}}{A tidy, long format data.frame of model
#'     predicted and observed CPUE by CPUE survey strata. This data.frame is now
#'     enhanced with new columns that include log-space standard error and
#'     associated confidence intervals that account for additional estimated
#'     observation error. If REMA is not run in multi-survey mode, or if CPUE
#'     data are not provided, an explanatory character string with instructions
#'     for fitting to CPUE data is returned.}
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
#' @seealso \code{\link{tidy_rema}}
#' @examples
#' \dontrun{
#' # placeholder for example
#' }
tidy_extra_cv <- function(tidy_rema,
                          save = FALSE,
                          path = NULL,
                          alpha_ci = 0.05) {
  # tidy_rema = tidy_rema(m1)
  # save = FALSE
  # path = NULL
  alpha_ci = 0.05


  if(nrow(tidy_rema$parameter_estimates %>% filter(parameter %in% c('extra_biomass_cv', 'extra_cpue_cv'))) == 0) {
    stop("Additional observation error for the biomass and/or CPUE index is not estimated in this model.")
  }

  if(nrow(tidy_rema$parameter_estimates %>% filter(parameter %in% c('extra_biomass_cv'))) > 1) {
    stop("The new tidy_extra_cv() and plot_extra_cv() functions only work when the additional estimated CV is shared across all biomass or CPUE strata.

         If you require this functionality, please file an issue at https://github.com/afsc-assessments/rema/issues")
  }

  if(nrow(tidy_rema$parameter_estimates %>% filter(parameter %in% c('extra_cpue_cv'))) > 1) {
    stop("The new tidy_extra_cv() and plot_extra_cv() functions only work when the additional estimated CV is shared across all biomass or CPUE strata.

         If you require this functionality, please file an issue at https://github.com/afsc-assessments/rema/issues")
  }

  if(nrow(tidy_rema$parameter_estimates %>% filter(parameter %in% c('extra_biomass_cv')) == 1)){

    tidy_rema$biomass_by_strata <- tidy_rema$biomass_by_strata %>%
      dplyr::mutate(extra_cv = tidy_rema$parameter_estimates %>%
                      filter(parameter == 'extra_biomass_cv') %>%
                      pull(estimate),
                    tot_sd_log_obs = ifelse(obs > 0, sqrt(log((obs_cv + extra_cv)^2 + 1)), NA),
                    tot_obs_lci = exp(log_obs - qnorm(1 - alpha_ci/2) * tot_sd_log_obs),
                    tot_obs_uci = exp(log_obs + qnorm(1 - alpha_ci/2) * tot_sd_log_obs))
  }

  if(nrow(tidy_rema$parameter_estimates %>% filter(parameter %in% c('extra_cpue_cv')) == 1)){

    tidy_rema$cpue_by_strata <- tidy_rema$cpue_by_strata %>%
      dplyr::mutate(extra_cv = tidy_rema$parameter_estimates %>%
                      filter(parameter == 'extra_cpue_cv') %>%
                      pull(estimate),
                    tot_sd_log_obs = ifelse(obs > 0, sqrt(log(obs_cv^2 + extra_cv^2 + 1)), NA),
                    tot_obs_lci = exp(log_obs - qnorm(1 - alpha_ci/2) * tot_sd_log_obs),
                    tot_obs_uci = exp(log_obs + qnorm(1 - alpha_ci/2) * tot_sd_log_obs))
  }

  return(tidy_rema)

}
