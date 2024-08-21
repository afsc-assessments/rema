#' Get one-step-head (OSA)
#'
#' Takes the rema model output from \code{\link{fit_rema}} and returns OSA
#' residuals calculated using
#' \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}} with accompanying
#' residual analysis plots. IMPORTANT: OSA residuals do not work for users
#' implementing the Tweedie distribution.
#'
#' @param rema_model list out output from \code{\link{fit_rema}}, which includes
#'   model results but also inputs. Of note to OSA residual calculations is the
#'   \code{rema_model$input$osa} object, which is a data.frame containing all
#'   the data or observations fit in the model that will have a residuals
#'   associated with them.
#' @param options list of options for calculating OSA residuals, passed to
#'   \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}. Default:
#'   \code{options = list(method = "fullGaussian", parallel = TRUE)}.
#'   Alternative methods include "cdf", "oneStepGeneric",
#'   "oneStepGaussianOffMode", and "oneStepGaussian".
#'
#' @return a list of tidied data.frames containing the biomass and CPUE survey
#'   residuals with accompanying data, as well as a QQ-plot, histogram
#'   of residuals, and plots of residuals~year and residuals~fitted values by
#'   strata for the biomass and CPUE survey.
#'
#' @import ggplot2
#' @import grid
#' @export
#' @seealso \code{\link{tidy_rema}}
#' @examples
#' \dontrun{
#' # placeholder for example
#' }
get_osa_residuals <- function(rema_model,
                              options = list(method = "fullGaussian",
                                             # "cdf","oneStepGeneric","oneStepGaussianOffMode","oneStepGaussian"),
                                             parallel = TRUE)) {

  # rema_model = m1
  # options = list(method = "fullGaussian", parallel = TRUE)

  print("**WARNING: OSA residuals do not currently work for the Tweedie distribution.**")

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

  p_qq <- function(dat, facet_strata = TRUE) {
    # dat = osa_resids
    sdnr <- sd(dat$residual)
    sdnr <- paste0('SDNR = ', formatC(sdnr, format = "f", digits = 2))
    p <- ggplot(data = dat, aes(sample = residual, col = strata)) +
      stat_qq() +
      # stat_qq_line() + # puts through IQR
      geom_abline(slope = 1, intercept = 0) +
      labs(x = 'Theoretical quantiles', y = 'Sample quantiles', col = 'Stratum') +
      ggplot2::scale_colour_viridis_d(direction = 1)

    if(isTRUE(facet_strata)) {p + facet_wrap(~strata)
    } else {p + annotation_custom(grid::grobTree(grid::textGrob(sdnr, x=0.1,  y=0.95, hjust=0,
                                              gp=grid::gpar(col="black", fontsize=10))))}

  }

  # p_hist <- function(dat) {
  #   ggplot(data = dat, aes(x = residual)) +
  #     geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white") +
  #     geom_density(fill = fill_alpha("#21918c", 0.6)) +
  #     labs(x = 'Residual', y = 'Density')
  # }

  p_fitted <- function(dat) {

    ybound <- max(abs(dat$residual), na.rm = TRUE) * 1.1

    ggplot(data = dat, aes(x = log_pred, y = residual)) +
      geom_hline(yintercept = 0, colour = "grey", size = 1) +
      geom_segment(aes(x = log_pred, xend = log_pred, y = 0, yend = residual)) +
      geom_point() +
      expand_limits(y = c(-ybound, ybound)) +
      facet_wrap(~strata, scales = 'free_x') +
      labs(x = 'Fitted values (log-scale)', y = 'Residual')
  }

  # p_acf <- function(dat) {
  #   # dat = osa_resids
  #   # dat = biomass_resids
  #
  #   biomacf <- acf(x = dat$residual, na.action = na.pass, plot = FALSE)
  #   acfci <- qnorm((1 + 0.95)/2)/sqrt(biomacf$n.used)
  #   biomacf <- data.frame(Lag = 1:(nrow(biomacf$acf)), ACF = biomacf$acf, lci = -acfci, uci = acfci)
  #   ggplot(biomacf, aes(x = Lag, y = ACF)) +
  #     geom_hline(yintercept = 0, colour = "grey", size = 1) +
  #     geom_hline(yintercept = biomacf$uci, colour = "blue", linetype = 2) +
  #     geom_hline(yintercept = biomacf$lci, colour = "blue", linetype = 2) +
  #     geom_segment(aes(x = Lag, xend = Lag, y = 0, yend = ACF)) +
  #     geom_point()
  #
  # }

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

      p1_qq <- p_qq(dat = osa_resids, facet_strata = FALSE)
      # p2_acf <- p_acf(dat = osa_resids)
      # p3_hist <- p_hist(dat = osa_resids)
      p1_biomass <- p_resids(dat = biomass_resids)
      p2_biomass <- p_fitted(dat = biomass_resids)
      p3_biomass <- p_qq(dat = biomass_resids, facet_strata = TRUE)

      biomass_resids <- biomass_resids %>%
        dplyr::select(model_name, variable, strata, year, log_pred, log_obs, residual)

      # if multi-survey
      if(is.data.frame(output$cpue_by_strata) & any(osa_resids$survey == 'CPUE survey')){

        cpue_resids <- output$cpue_by_strata %>%
          dplyr::left_join(osa_resids %>%
                             dplyr::filter(survey == 'CPUE survey') %>%
                             dplyr::select(year, strata, residual))

        p1_cpue <- p_resids(dat = cpue_resids)
        p2_cpue <- p_fitted(dat = cpue_resids)
        p3_cpue <- p_qq(dat = cpue_resids, facet_strata = TRUE)

        cpue_resids <- cpue_resids %>%
          dplyr::select(model_name, variable, strata, year, log_pred, log_obs, residual)
      }

      out <- list(residuals = NULL,
                  plots = NULL)

      out$residuals <- list(biomass = biomass_resids,
                            cpue = cpue_resids)

      out$plots <- list(qq = p1_qq,
                        # acf = p2_acf,
                        # histo = p3_hist,
                        biomass_resids = p1_biomass,
                        biomass_fitted = p2_biomass,
                        biomass_qq = p3_biomass,
                        cpue_resids = p1_cpue,
                        cpue_fitted = p2_cpue,
                        cpue_qq = p3_cpue)
      return(out)
      }

  } else warning(paste("","** Did not do OSA residual analyses. **",
                       "Error during TMB::sdreport(). Check for unidentifiable parameters.","",sep='\n'))

}
