#' Fit REMA model
#'
#' Fits the compiled REMA model using
#' \code{\link[TMB:MakeADFun]{TMB::MakeADFun}} and
#' \code{\link[stats:nlminb]{stats::nlminb}}. Source code and documentation modified from
#' \href{https://github.com/timjmiller/wham/blob/master/R/fit_wham.R}{wham::fit_wham}.
#'
#' Future development: Implement one-step-ahead (OSA) residuals for evaluating
#' model goodness-of-fit \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}).
#' OSA residuals are more appropriate than standard residuals for models with
#' random effects
#' (\href{https://link.springer.com/article/10.1007/s10651-017-0372-4}{Thygeson
#' et al. (2017)}. See
#' \href{https://github.com/timjmiller/wham/blob/master/R/fit_wham.R}{wham} for
#' an example of OSA implementation and additional OSA residual options (e.g.
#' full Gaussian approximation instead of the (default) generic method using
#' \code{osa.opts=list(method="fullGaussian")}.
#'
#' @param input Named list output from \code{\link{prepare_rema_input}}, which
#'   includes the following components needed to fit model using
#'   \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}:
#'   \describe{
#'     \item{\code{$data}}{Data, a list of data objects for model fitting or
#'     specification (e.g., user-defined pentalties, index pointers, etc.). A
#'     required input to \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{$par}}{Parameters, a list of all random and fixed effects
#'     parameter objects.  A required input to \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{$map}}{Map, a mechanism for collecting and fixing parameters
#'     in TMB.  An input to \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{$random}}{Character vector defining the parameters to treat
#'     as random effects. An input to \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{$model_name}}{Character, name of the model, e.g. \code{"GOA
#'     shortraker with LLS by depth strata"}. Useful for model comparison.}
#'   }
#' @param n.newton integer, number of additional Newton steps after
#'   optimization. Not an option that is currently needed, but is passed to
#'   \code{\link{fit_tmb}}. Default = \code{0}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See
#'   \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#' @param model (optional), a previously fit rema model.
#' @param do.check T/F, check if model parameters are identifiable? Passed to
#'   \code{\link{fit_tmb}}. Runs internal function \code{check_estimability},
#'   originally provided by https://github.com/kaskr/TMB_contrib_R/TMBhelper.
#'   Default = \code{TRUE}.
#' @param MakeADFun.silent T/F, Passed to silent argument of
#'   \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{TRUE}.
#' @param do.fit T/F, fit the model using \code{fit_tmb}. Default = \code{TRUE}.
#' @param save.sdrep T/F, save the full \code{\link[TMB]{TMB::sdreport}} object?
#'   If \code{FALSE}, only save
#'   \code{\link[TMB:summary.sdreport]{summary.sdreport}} to reduce model object
#'   file size. Default = \code{TRUE}.
#'
#' @return a fit TMB model with additional output if specified:
#'   \describe{
#'     \item{\code{$rep}}{List of derived quantity estimates (e.g. estimated
#'     biomass)}
#'     \item{\code{$sdrep}}{Parameter estimates (and standard errors if
#'     \code{do.sdrep = TRUE})}
#'   }
#'
#' @useDynLib rema
#' @export
#'
#' @seealso \code{\link{fit_tmb}},
#'          \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}}
#'
#' @examples
#' \dontrun{
#' # place holder for example code
#' }
fit_rema <- function(input,
                     n.newton = 0,
                     do.sdrep = TRUE,
                     model = NULL,
                     do.check = FALSE,
                     MakeADFun.silent = TRUE,
                     do.fit = TRUE,
                     save.sdrep = TRUE) {
  # n.newton = 1
  # do.sdrep = TRUE
  # do.osa = TRUE
  # osa.opts = list(method = "cdf", parallel = TRUE)
  # model = NULL
  # do.check = FALSE
  # MakeADFun.silent = TRUE
  # do.fit = TRUE
  # save.sdrep = TRUE

  # fit model
  if(is.null(model)){
    mod <- TMB::MakeADFun(input$data, input$par, DLL = "rema", random = input$random, map = input$map, silent = MakeADFun.silent)
  } else {mod = model}

  mod$years <- input$data$model_yrs
  mod$model_name <- input$model_name
  mod$input <- input
  ver <- sessioninfo::package_info() %>% as.data.frame %>% dplyr::filter(package == "rema") %>% dplyr::select(loadedversion, source) %>% unname
  mod$rema_version <- paste0(ver, collapse=" / ")
  if(do.fit){
    btime <- Sys.time()
    mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = FALSE, do.check = do.check, save.sdrep = save.sdrep)
    mod$runtime <- round(difftime(Sys.time(), btime, units = "s"), 1)

    # SEs for estimated parameters. Variance calculations use the delta method.
    if(do.sdrep) {
      mod$sdrep <- try(TMB::sdreport(mod))
      mod$is_sdrep <- !is.character(mod$sdrep)
      if(mod$is_sdrep) mod$na_sdrep <- any(is.na(summary(mod$sdrep,"fixed")[,2])) else mod$na_sdrep = NA
      if(!save.sdrep) mod$sdrep <- summary(mod$sdrep) # only save summary to reduce model object size
      check_convergence(mod)
    }

    # error message reporting
    if(!is.null(mod$err)) warning(paste("","** Error during model fit. **",
                                        "Check for unidentifiable parameters.","",mod$err,"",sep='\n'))
  }
  else { # model not fit, but generate report and parList without fitted model. potential use for future development (e.g. projections?)
    mod$rep = mod$report() # par values don't matter because function has not been evaluated
    mod$parList = mod$env$parList()
    warning("Model not fit. Report and parList not valid because objective function has not been evaluated. Try do.fit = TRUE.")
  }

  return(mod)
}
