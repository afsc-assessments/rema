#' Plot survey data and model output
#'
#' Takes list output from \code{\link{tidy_rema}} and returns a list of
#' \code{ggplot2} objects to be plotted or saved.
#'
#' @param tidy_rema list out output from \code{\link{tidy_rema}}, which includes
#'   model results but also inputs
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
#' @import ggplot2
#' @export
#' @seealso \code{\link{tidy_rema}}
#' @examples
#' \dontrun{
#' # placeholder for example
#' }
plot_rema <- function(tidy_rema,
                      save = FALSE,
                      filetype = "png",
                      path = NULL,
                      xlab = NULL,
                      biomass_ylab = 'Biomass',
                      cpue_ylab = 'CPUE') {
  # DELETE ME
  # tidy_rema = m_output
  # save = FALSE
  # filetype = "png"
  # path = NULL
  # xlab = NULL
  # biomass_ylab = 'Biomass (t)'
  # cpue_ylab = 'Relative Population Weight'

  rema_plots <- list()

  # how many plot rows should I have?

  if(!is.data.frame(tidy_rema$biomass_by_strata)) {
    stop("Something went wrong... Please review the ?prepare_rema_input, ?fit_rema, and make sure REMA converged by running check_convergence(my_rema_model).")
  } else {
    p1 <- ggplot(data = tidy_rema$biomass_by_strata,
                 aes(x = year, y = pred)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                  col = 'grey', fill = 'grey') +
      geom_line() +
      facet_wrap(~strata, nrow = NULL) +
      geom_point(aes(x = year, y = obs)) +
      geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci)) +
      scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = biomass_ylab)
  }

  if(!is.data.frame(tidy_rema$cpue_by_strata)) {
    p2 <- "REMA was fit only to biomass survey data, therefore no CPUE predictions available. If the user has a CPUE survey index and wants to fit to it, please see ?prepare_rema_input() for details."
  } else {
    p2 <- ggplot(data = tidy_rema$cpue_by_strata,
                 aes(x = year, y = pred)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                  col = 'grey', fill = 'grey') +
      geom_line() +
      facet_wrap(~strata, nrow = NULL) +
      geom_point(aes(x = year, y = obs)) +
      geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci)) +
      scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = cpue_ylab)
  }

  if(!is.data.frame(tidy_rema$biomass_by_cpue_strata)) {
    p3 <- "'biomass_by_cpue_strata' is reserved for multi-survey scenarios when there are more biomass survey strata than CPUE survey strata, and the user wants predicted biomass at the same resolution as the CPUE survey index."
  } else {
    p3 <- ggplot(data = tidy_rema$biomass_by_cpue_strata,
                 aes(x = year, y = pred)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                  col = 'grey', fill = 'grey') +
      geom_line() +
      facet_wrap(~strata, nrow = NULL) +
      # geom_point(aes(x = year, y = obs)) +
      # geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci)) +
      scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = biomass_ylab)
  }

  if(!is.data.frame(tidy_rema$total_predicted_biomass)) {
    stop("Something went wrong... Please review the ?prepare_rema_input, ?fit_rema, and make sure REMA converged by running check_convergence(my_rema_model).")
  } else {

    p4 <- ggplot(data = tidy_rema$total_predicted_biomass,
                 aes(x = year, y = pred)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                  col = 'grey', fill = 'grey') +
      geom_line() +
      scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = biomass_ylab)
  }

  if(!is.data.frame(tidy_rema$total_predicted_cpue)) {
    p5 <- "REMA was fit only to biomass survey data, therefore no CPUE predictions available. If the user has a CPUE survey index and wants to fit to it, please see ?prepare_rema_input() for details."
  } else {

    p5 <- ggplot(data = tidy_rema$total_predicted_cpue,
                 aes(x = year, y = pred)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                  col = 'grey', fill = 'grey') +
      geom_line() +
      scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = cpue_ylab)
  }

  rema_plots$biomass_by_strata <- p1
  rema_plots$cpue_by_strata <- p2
  rema_plots$biomass_by_cpue_strata <- p3
  rema_plots$total_predicted_biomass <- p4
  rema_plots$total_predicted_cpue <- p5

  return(rema_plots)

}
