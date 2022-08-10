#' Plot REMA model comparisons and return AIC values when appropriate
#'
#' Takes list of REMA models from from \code{\link{fit_rema}}, and returns a
#' list of \code{ggplot2} objects to be plotted or saved, a list of
#' \code{\link{tidy_rema}} data.frames, and AIC values.
#'
#' @param rema_models list of REMA models to be compared. Each REMA model in the
#'   list should be a list object output from \code{\link{fit_rema}}
#' @param admb_re list of ADMB RE model input/output from
#'   \code{\link{read_admb_re}}. Accepts a single list, not list of multiple
#'   ADMB RE models. If admb_re is provided, no AIC calculations will be
#'   conducted.
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
#' @return a list with the following items:
#'   \describe{
#'     \item{\code{$output}}{A list of tidied dataframes that include parameter
#'     estimates, biomass and optional CPUE data, and REMA model predictions for
#'     each model to be compared. Results for a given variable are only included
#'     if they are applicable to all comparison models. For example, if CPUE is
#'     fit in one model but not another, \code{compare$output$cpue_by_strata)}
#'     will return an informational message instead of a dataframe. See
#'     \code{\link{tidy_rema}} for more information.}
#'     \item{\code{$plots}}{ggplot2 figure objects of \code{compare$output} data.}
#'     \item{\code{$aic}}{A dataframe of Akaike Information Criteria (AIC)
#'     values. Only output if the underlying models are fit to the same data.}
#' }
#'
#' @import ggplot2
#' @export
#' @seealso \code{\link{tidy_rema}}, \code{\link{plot_rema}}
#' @examples
#' \dontrun{
#' # placeholder for example
#' }
compare_rema_models <- function(rema_models,
                                admb_re = NULL,
                                save = FALSE,
                                filetype = "png",
                                path = NULL,
                                xlab = NULL,
                                biomass_ylab = 'Biomass',
                                cpue_ylab = 'CPUE') {
  # rema_models <- list(m, m2)
  # admb_re = NULL
  # biomass_ylab <- 'ROV biomass'
  # cpue_ylab <- 'IPHC setline survey CPUE'
  # xlab = NULL

  if(!is.list(rema_models)) {
    stop("This function expects a list of model objects (e.g., compare_rema_models(rema_models = list(m1, m2))). Each of the model objects should be a list returned from fit_rema() See ?compare_rema_models for details.")
  }

  compare_rema_output <- list()
  compare_rema_plots <- list()
  compare_rema_aic <- list()

  # tidy out from all the selected models
  out <- lapply(rema_models, tidy_rema)

  if(!is.null(admb_re)) {
    warning("The biomass data in admb_re_results$biomass_by_strata assumes that any zero biomass observations were removed in the model fitting process (i.e. assumed to be NAs or failed surveys), which is the most common assumption across RE.tpls. Please check this assumption in the RE.tpl before using admb_re$admb_re_results to compare REMA model fits to the data.")
    out[[length(out)+1]] <- admb_re$admb_re_results
  }

  # pull unique variables from each model output list into its own list
  parameter_estimates <- lapply(out, `[[`, 1)
  biomass_by_strata <- lapply(out, `[[`, 2)
  cpue_by_strata <- lapply(out, `[[`, 3)
  biomass_by_cpue_strata <- lapply(out, `[[`, 4)
  total_predicted_biomass <- lapply(out, `[[`, 5)
  total_predicted_cpue <- lapply(out, `[[`, 6)

  # check if each model output data is structured the same or if some have
  # messages as output (e.g. one model has CPUE and the other doesn't, so one
  # model will output a data.frame for cpue_by_strata but the other will output
  # a character string msg)
  tst_parameter_estimates <- all(unlist(lapply(parameter_estimates, is.data.frame)))
  tst_biomass_by_strata <- all(unlist(lapply(biomass_by_strata, is.data.frame)))
  tst_cpue_by_strata <- all(unlist(lapply(cpue_by_strata, is.data.frame)))
  tst_biomass_by_cpue_strata <- all(unlist(lapply(biomass_by_cpue_strata, is.data.frame)))
  tst_total_predicted_biomass <- all(unlist(lapply(total_predicted_biomass, is.data.frame)))
  tst_total_predicted_cpue <- all(unlist(lapply(total_predicted_cpue, is.data.frame)))

  # output tidy parameter estimates
  if(tst_parameter_estimates) {
    out_parameter_estimates <- do.call('rbind', parameter_estimates)
  } else if(!is.null(admb_re)){
    out_parameter_estimates <- "Parameter estimates for the ADMB version of the RE model are not readily available for comparison with REMA models. User will need to check the ADMB .std file for parameter estimates."
  } else {
    stop("Something went wrong... run check_convergence() and review output$parameter_estimates from tidy_rema() output for all REMA models you want to compare. All models should have should meet minimum convergence criteria and have valid parameter estimates.")
  }

  # biomass by strata
  if(tst_biomass_by_strata) {

    biomass_by_strata <- lapply(biomass_by_strata, function(x) {x %>%
        dplyr::select(model_name, strata, variable, year, pred, pred_lci, pred_uci, obs, obs_cv, obs_lci, obs_uci)})

    out_biomass_by_strata <- do.call('rbind', biomass_by_strata)

    # test that biomass data are all equal
    tst_biomass_data <- lapply(biomass_by_strata, `[`, 8:9) # 8:9 = obs and obs_cv
    tst_biomass_data <- all(sapply(tst_biomass_data, identical, tst_biomass_data[[1]]))
    tst_biomass_data <- TRUE

    if(isFALSE(tst_biomass_data)) {
      out_biomass_by_strata <- "The REMA models selected for comparison were fit to different biomass data, and therefore the fits to the biomass data by strata cannot be compared."
      p1 <- "The REMA models selected for comparison were fit to different biomass data, and therefore the fits to the biomass data by strata cannot be compared."
    }

    if(isTRUE(tst_biomass_data)) {

      p1 <- ggplot(data = out_biomass_by_strata,
                   aes(x = year, y = pred,
                       col = model_name)) +
        geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                        fill = model_name), col = NA,
                    alpha = 0.25) +
        geom_line() +
        facet_wrap(~strata, nrow = NULL) +
        # geom_point(aes(x = year, y = obs), col = 'black') +
        # geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci), col = 'black') +
        geom_point(aes(x = year, y = obs, col = model_name, shape = model_name)) +
        geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci, col = model_name)) +
        scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
        labs(x = xlab, y = biomass_ylab,
             fill = NULL, colour = NULL, shape = NULL) +
        # ggplot2::scale_colour_brewer(palette = 'Set1') +
        # ggplot2::scale_fill_brewer(palette = 'Set1')
        ggplot2::scale_fill_viridis_d(direction = 1) +
        ggplot2::scale_colour_viridis_d(direction = 1)
    }
  } else if(!is.null(admb_re)) {
    out_biomass_by_strata <- "Biomass estimates for the ADMB version of the RE model do not appear to be readily available for comparison with REMA models. Check the rwout.rep file and ?read_admb_re for more information."
    p1 <- "Biomass estimates for the ADMB version of the RE model do not appear to be readily available for comparison with REMA models. Check the rwout.rep file and ?read_admb_re for more information."
  } else {
    stop("Something went wrong... run check_convergence() and review output$biomass_by_strata from tidy_rema() output for all REMA models you want to compare. All models should have should meet minimum convergence criteria and have valid biomass data and model predictions.")
  }

  # cpue by strata
  if(tst_cpue_by_strata) {

    cpue_by_strata <- lapply(cpue_by_strata, function(x) {x %>%
        dplyr::select(model_name, strata, variable, year, pred, pred_lci, pred_uci, obs, obs_cv, obs_lci, obs_uci)})

    out_cpue_by_strata <- do.call('rbind', cpue_by_strata)

    # test that cpue data are all equal
    tst_cpue_data <- lapply(cpue_by_strata, `[`, 8:9) # 8:9 = obs and obs_cv
    tst_cpue_data <- all(sapply(tst_cpue_data, identical, tst_cpue_data[[1]]))
    tst_cpue_data <- TRUE

    if(isFALSE(tst_cpue_data)) {
      out_cpue_by_strata <- "The REMA models selected for comparison were fit to different CPUE data, and therefore the fits to the CPUE data by strata cannot be compared."
    }
    if(isTRUE(tst_cpue_data)) {

      p2 <- ggplot(data = out_cpue_by_strata,
                   aes(x = year, y = pred,
                       col = model_name)) +
        geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                        fill = model_name), col = NA,
                    alpha = 0.25) +
        geom_line() +
        facet_wrap(~strata, nrow = NULL) +
        # geom_point(aes(x = year, y = obs), col = 'black') +
        # geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci), col = 'black') +
        geom_point(aes(x = year, y = obs, col = model_name, shape = model_name)) +
        geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci, col = model_name)) +
        scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
        labs(x = xlab, y = cpue_ylab,
             fill = NULL, colour = NULL, shape = NULL) +
        ggplot2::scale_fill_viridis_d(direction = 1) +
        ggplot2::scale_colour_viridis_d(direction = 1) #+
        # ggplot2::scale_colour_brewer(palette = 'Set1') +
        # ggplot2::scale_fill_brewer(palette = 'Set1')
    }
  } else {
    out_cpue_by_strata <- "One or more of the models selected for comparison were not fit to CPUE data and therefore cannot be compared."
    p2 <- "One or more of the models selected for comparison were not fit to CPUE data and therefore cannot be compared."
  }

  # biomass by cpue strata
  if(tst_biomass_by_cpue_strata) {
    out_biomass_by_cpue_strata <- do.call('rbind', biomass_by_cpue_strata)

    p3 <- ggplot(data = out_biomass_by_cpue_strata,
                 aes(x = year, y = pred,
                     col = model_name)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                      fill = model_name), col = NA,
                  alpha = 0.25) +
      geom_line() +
      facet_wrap(~strata, nrow = NULL) +
      # geom_point(aes(x = year, y = obs), col = 'black') +
      # geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci), col = 'black') +
      geom_point(aes(x = year, y = obs, col = model_name, shape = model_name)) +
      geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci, col = model_name)) +
      scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = biomass_ylab,
           fill = NULL, colour = NULL, shape = NULL) +
      ggplot2::scale_fill_viridis_d(direction = 1) +
      ggplot2::scale_colour_viridis_d(direction = 1)
      # ggplot2::scale_colour_brewer(palette = 'Set1') +
      # ggplot2::scale_fill_brewer(palette = 'Set1')

  } else {
    out_biomass_by_cpue_strata <- "'biomass_by_cpue_strata' is reserved for multi-survey scenarios when there are more biomass survey strata than CPUE survey strata, and the user wants predicted biomass at the same resolution as the CPUE survey index. One or more of the models selected for comparison did not meet this criterion."
    p3 <- "'biomass_by_cpue_strata' is reserved for multi-survey scenarios when there are more biomass survey strata than CPUE survey strata, and the user wants predicted biomass at the same resolution as the CPUE survey index. One or more of the models selected for comparison did not meet this criterion."
  }

  # total predicted biomass
  if(tst_total_predicted_biomass) {

    total_predicted_biomass <- lapply(total_predicted_biomass, function(x) {x %>%
        dplyr::select(model_name, variable, year, pred, pred_lci, pred_uci)})

    out_total_predicted_biomass <- do.call('rbind', total_predicted_biomass)

    p4 <- ggplot(data = out_total_predicted_biomass,
                 aes(x = year, y = pred,
                     col = model_name)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                      fill = model_name), col = NA,
                  alpha = 0.25) +
      geom_line() +
      scale_y_continuous(labels = scales::comma) + #, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = biomass_ylab,
           fill = NULL, colour = NULL) +
      ggplot2::scale_fill_viridis_d(direction = 1) +
      ggplot2::scale_colour_viridis_d(direction = 1)
      # ggplot2::scale_colour_brewer(palette = 'Set1') +
      # ggplot2::scale_fill_brewer(palette = 'Set1')
  } else if(!is.null(admb_re)) {
    out_total_predicted_biomass <- "Biomass estimates for the ADMB version of the RE model do not appear to be readily available for comparison with REMA models. Check the rwout.rep file and ?read_admb_re for more information."
    p4 <- "Biomass estimates for the ADMB version of the RE model do not appear to be readily available for comparison with REMA models. Check the rwout.rep file and ?read_admb_re for more information."
  } else {
    stop("Something went wrong... run check_convergence() and review output$total_predicted_biomass from tidy_rema() output for all REMA models you want to compare. All models should have should meet minimum convergence criteria and have valid model predictions of total biomass.")
  }

  # total predicted cpue
  if(tst_total_predicted_cpue) {

    total_predicted_cpue <- lapply(total_predicted_cpue, function(x) {x %>%
        dplyr::select(model_name, variable, year, pred, pred_lci, pred_uci)})

    out_total_predicted_cpue <- do.call('rbind', total_predicted_cpue)

    p5 <- ggplot(data = out_total_predicted_cpue,
                 aes(x = year, y = pred,
                     col = model_name)) +
      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                      fill = model_name), col = NA,
                  alpha = 0.25) +
      geom_line() +
      scale_y_continuous(labels = scales::comma) + #, expand = c(0, 0), limits = c(0, NA)) +
      labs(x = xlab, y = cpue_ylab,
           fill = NULL, colour = NULL) +
      ggplot2::scale_fill_viridis_d(direction = 1) +
      ggplot2::scale_colour_viridis_d(direction = 1)
      # ggplot2::scale_colour_brewer(palette = 'Set1') +
      # ggplot2::scale_fill_brewer(palette = 'Set1')

  } else {
    out_total_predicted_cpue <- "Either one or more of the models selected for comparison were not fit to CPUE survey data OR the CPUE survey index was defined as not summable in prepare_rema_input(). If the CPUE index is summable (e.g. Relative Population Numbers), please select sum_cpue_index = TRUE in prepare_rema_input(). See ?prepare_rema_input() for more details."
    p5 <- "Either one or more of the models selected for comparison were not fit to CPUE survey data OR the CPUE survey index was defined as not summable in prepare_rema_input(). If the CPUE index is summable (e.g. Relative Population Numbers), please select sum_cpue_index = TRUE in prepare_rema_input(). See ?prepare_rema_input() for more details."
  }

  # AIC calculations
  out_aic <- NULL
  run_aic <- NULL

  # no AIC for admb models
  if(!is.null(admb_re)) {
    run_aic <- FALSE
    out_aic <- 'AIC calculations are not currently available for ADMB RE models.'

    # all models to be compared are run with multi_survey = 1, check that the
    # biomass and survey data are equal before calculating AIC
  } else if(all(sapply(rema_models, function(x) {isTRUE(x$input$data$multi_survey == 1)})) &
            isTRUE(tst_biomass_data) & !is.vector(p2)) # p2 = tst_cpue_data
   {
    run_aic <- TRUE
    # if all models are run with multi_survey = 0, check that biomass data are
    # equal before calculating AIC
  } else if(all(sapply(rema_models, function(x) {isTRUE(x$input$data$multi_survey == 0)})) &
            isTRUE(tst_biomass_data)
  ) {
    run_aic <- TRUE
    # otherwise do not calculate aic
  } else {
    run_aic <- FALSE
    out_aic <- 'AIC could not be calculated because the REMA models selected for comparison were fit to different biomass or CPUE survey data, have different assumptions about observed zero biomass or CPUE data, or some models were fit in multi-survey model while others were not.'
  }

  # run AIC calcs when appropriate
  if(isTRUE(run_aic)) {

    aic <- sapply(rema_models, function(x) {
      k = length(x$opt$par)
      2*(x$opt$obj + k)
    })
    model_names <- sapply(rema_models, function(x) {x$input$model_name})
    aic <- round(aic, 1)
    daic <- round(aic - min(aic), 1)

    out_aic <- data.frame(model_name = model_names,
                          aic = aic,
                          daic = daic) %>%
      dplyr::arrange(daic)

  }

  # Prepare function output
  compare_rema_output$parameter_estimates <- out_parameter_estimates
  compare_rema_output$biomass_by_strata <- out_biomass_by_strata
  compare_rema_output$cpue_by_strata <- out_cpue_by_strata
  compare_rema_output$biomass_by_cpue_strata <- out_biomass_by_cpue_strata
  compare_rema_output$total_predicted_biomass <- out_total_predicted_biomass
  compare_rema_output$total_predicted_cpue <- out_total_predicted_cpue

  compare_rema_plots$biomass_by_strata <- p1
  compare_rema_plots$cpue_by_strata <- p2
  compare_rema_plots$biomass_by_cpue_strata <- p3
  compare_rema_plots$total_predicted_biomass <- p4
  compare_rema_plots$total_predicted_cpue <- p5

  compare_rema_aic <- out_aic

  fout <- list(output = compare_rema_output,
               plots = compare_rema_plots,
               aic = compare_rema_aic)
  return(fout)
}

