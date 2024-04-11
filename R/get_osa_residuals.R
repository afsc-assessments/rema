#' Get one-step-head (OSA)
#'
#' Takes list output from \code{\link{tidy_rema}} and returns a list of
#' \code{ggplot2} objects to be plotted or saved. This feature is
#' experimental, and OSA residuals are still in Beta mode in
#' \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}. The default "cdf"
#' method sometimes results in NA values for residuals, especially when
#' observation errors are small.
#'
#' @param rema_model list out output from \code{\link{fit_rema}}, which includes
#'   model results but also inputs. Of note to OSA residual calculations is the
#'   \code{rema_model$input$osa} object, which is a data.frame containing all
#'   the data or observations fit in the model that will have a residuals
#'   associated with them.
#' @param options list of options for calculating OSA residuals, passed to
#'   \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}. Default:
#'   \code{options = list(method = "cdf", parallel = TRUE)}.
#'
#' @return a list of tidied data.frames containing the biomass and CPUE survey
#'   residuals with accompanying data, as well as ggplot2 plots to visual
#'   results from the OSA residual analysis.
#'
#' @import ggplot2
#' @export
#' @seealso \code{\link{tidy_rema}}
#' @examples
#' \dontrun{
#' # placeholder for example
#' }
get_osa_residuals <- function(rema_model,
                              options = list(method =  "oneStepGeneric",
                                             # "cdf","oneStepGaussianOffMode",
                                             # , "fullGaussian",
                                             # "oneStepGaussian"),
                                             parallel = TRUE)) {

  # rema_model = m
  # options = list(method = "cdf", parallel = TRUE)


  p_resids <- function(dat) {

    ybound <- max(abs(dat$residual), na.rm = TRUE) * 1.1

    ggplot(data = dat, aes(x = year, y = residual)) +
      geom_hline(yintercept = 0, colour = "grey", size = 1) +
      geom_segment(aes(x = year, xend = year, y = 0, yend = residual)) +
      geom_point() +
      expand_limits(y = c(-ybound, ybound)) +
      facet_wrap(~strata) +
      labs(x = 'Year', y = 'Residual')
  }

  p_qq <- function(dat) {
    ggplot(data = dat, aes(sample = residual)) +
      stat_qq() +
      stat_qq_line() +
      facet_wrap(~strata) +
      labs(x = 'Theoretical quantiles', y = 'Sample quantiles')
  }

  p_hist <- function(dat) {
    ggplot(data = dat, aes(x = residual)) +
      # geom_histogram() +
      geom_histogram(aes(y = ..density..), colour = "black", fill = "white")+
      geom_density(alpha = 0.2, fill = "#FF6666") +
      facet_wrap(~strata, scales = 'free') +
      labs(x = 'Residual', y = 'Density')
  }

  p_fitted <- function(dat) {

    ybound <- max(abs(dat$residual), na.rm = TRUE) * 1.1

    ggplot(data = dat, aes(x = log_pred, y = residual)) +
      geom_hline(yintercept = 0, colour = "grey", size = 1) +
      geom_point() +
      expand_limits(y = c(-ybound, ybound)) +
      facet_wrap(~strata, scales = 'free_x') +
      labs(x = 'Fitted values (log-scale)', y = 'Residual')
  }

  if(is.null(rema_model$err)) {

    if(rema_model$is_sdrep) { # only do OSA residuals if sdrep ran
      cat("Doing OSA residuals...\n")

      osa_resids <- rema_model$input$osa
      osa_resids$residual <- NA
      r <- suppressWarnings(TMB::oneStepPredict(obj = rema_model,
                                                observation.name = "obsvec",
                                                data.term.indicator = "keep",
                                                discrete = FALSE,
                                                method = options$method))
      osa_resids$residual <- r$residual

      # link residuals to observed and predicted values
      output <- tidy_rema(rema_model)

      biomass_resids <- p1_biomass <-  p2_biomass <-  p3_biomass <- p4_biomass <- NULL
      cpue_resids <- p1_cpue <- p2_cpue <- p3_cpue <- p4_cpue <- NULL

      biomass_resids <- output$biomass_by_strata %>%
        dplyr::left_join(osa_resids %>%
                           dplyr::filter(survey == 'Biomass survey') %>%
                           dplyr::select(year, strata, residual))

      p1_biomass <- p_resids(dat = biomass_resids)
      p2_biomass <- p_qq(dat = biomass_resids)
      p3_biomass <- p_hist(dat = biomass_resids)
      p4_biomass <- p_fitted(dat = biomass_resids)

      biomass_resids <- biomass_resids %>%
        dplyr::select(model_name, variable, strata, year, log_pred, log_obs, residual)

      # if multi-survey
      if(is.data.frame(output$cpue_by_strata) & any(osa_resids$survey == 'CPUE survey')){

        cpue_resids <- output$cpue_by_strata %>%
          dplyr::left_join(osa_resids %>%
                             dplyr::filter(survey == 'CPUE survey') %>%
                             dplyr::select(year, strata, residual))

        p1_cpue <- p_resids(dat = cpue_resids)
        p2_cpue <- p_qq(dat = cpue_resids)
        p3_cpue <- p_hist(dat = cpue_resids)
        p4_cpue <- p_fitted(dat = cpue_resids)

        cpue_resids <- cpue_resids %>%
          dplyr::select(model_name, variable, strata, year, log_pred, log_obs, residual)
      }

      out <- list(residuals = NULL,
                  plots = NULL)

      out$residuals <- list(biomass = biomass_resids,
                            cpue = cpue_resids)

      out$plots <- list(biomass_resids = p1_biomass,
                        biomass_qqplot = p2_biomass,
                        biomass_hist = p3_biomass,
                        biomass_fitted = p4_biomass,
                        cpue_resids = p1_cpue,
                        cpue_qqplot = p2_cpue,
                        cpue_hist = p3_cpue,
                        cpue_fitted = p4_cpue)
      return(out)
      }

  } else warning(paste("","** Did not do OSA residual analyses. **",
                       "Error during TMB::sdreport(). Check for unidentifiable parameters.","",sep='\n'))

}
