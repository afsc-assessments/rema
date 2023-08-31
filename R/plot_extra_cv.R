#' Plot the additional estimated observation error for biomass by strata and/or cpue by strata
#'
#' Takes list output from \code{\link{tidy_rema}} and returns a list of
#' \code{ggplot2} objects to be plotted or saved.
#'
#' @param tidy_rema list out output from \code{\link{tidy_extra_cv}}, which
#'   includes inputs, model results, and confidence intervals for the total
#'   observation error (fixed + estimated)
#' @param save (optional) logical (T/F) save figures as \code{filetype} in
#'   \code{path}. Default = FALSE. NOT YET IMPLEMENTED.
#' @param filetype (optional) character string; type of figure file. Default =
#'   'png'. NOT YET IMPLEMENTED.
#' @param path (optional) directory path to location where figure files are to
#'   be saved if \code{save = TRUE}. NOT YET IMPLEMENTED.
#' @param xlab (optional) label for x-axis of biomass and CPUE plots (e.g.
#'   'Year'). Default = NULL.
#' @param biomass_ylab (optional) label for y-axis of biomass plots (e.g.
#'   'Biomass (t)'). Default = 'Biomass'.
#' @param cpue_ylab (optional) label for y-axis of CPUE plots (e.g. 'Relative
#'   Population Number'). Default = 'CPUE'.
#'
#' @return a list of ggplot2 plots or character string messages about the data.
#'   Except for parameter estimates, the objects output from
#'   \code{\link{tidy_rema}} are the same outputted from this function.
#'
#' @import ggplot2
#' @export
#' @seealso \code{\link{tidy_rema}}
#' @examples
#' \dontrun{
#' # placeholder for example
#' }
plot_extra_cv <- function(tidy_rema,
                          save = FALSE,
                          filetype = "png",
                          path = NULL,
                          xlab = NULL,
                          biomass_ylab = 'Biomass',
                          cpue_ylab = 'CPUE') {
  # tidy_rema = tidy_rema(m1)
  # save = FALSE
  # filetype = "png"
  # path = NULL
  # xlab = NULL
  # biomass_ylab = 'Biomass (t)'
  # cpue_ylab = 'Relative Population Weight'

  rema_plots <- list()

  if(!is.data.frame(tidy_rema$biomass_by_strata)) {
    stop("Something went wrong... Please review the ?prepare_rema_input, ?fit_rema, and make sure REMA converged by running check_convergence(my_rema_model).")
  } else if("tot_obs_lci" %in% colnames(tidy_rema$biomass_by_strata)) {
    p1 <- ggplot(data = tidy_rema$biomass_by_strata,
                 aes(x = year, y = pred)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                  col = 'grey', fill = 'grey') +
      geom_line() +
      facet_wrap(~strata, nrow = NULL) +
      geom_point(aes(x = year, y = obs)) +
      geom_errorbar(aes(x = year, ymin = tot_obs_lci, ymax = tot_obs_uci), position = position_dodge(0), width = 0) +
      geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci), linewidth = 1) +
      scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = biomass_ylab)
  } else {
    p1 <- "Additional observation error is not estimated for biomass in this model. Please use plot_rema() instead."
  }

  if(!is.data.frame(tidy_rema$cpue_by_strata)) {
    p2 <- "REMA was fit only to biomass survey data, therefore no CPUE predictions available. If the user has a CPUE survey index and wants to fit to it, please see ?prepare_rema_input() for details."
  } else if("tot_obs_lci" %in% colnames(tidy_rema$cpue_by_strata)) {
    p2 <- ggplot(data = tidy_rema$cpue_by_strata,
                 aes(x = year, y = pred)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                  col = 'grey', fill = 'grey') +
      geom_line() +
      facet_wrap(~strata, nrow = NULL) +
      geom_point(aes(x = year, y = obs)) +
      geom_errorbar(aes(x = year, ymin = tot_obs_lci, ymax = tot_obs_uci), position = position_dodge(0), width = 0) +
      geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci), linewidth = 1) +
      scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = cpue_ylab)
  } else {
    p2 <- "Additional observation error is not estimated for CPUE in this model. Please use plot_rema() instead."
  }

  rema_plots$biomass_by_strata <- p1
  rema_plots$cpue_by_strata <- p2

  return(rema_plots)

}
