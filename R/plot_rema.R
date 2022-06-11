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
  tidy_rema =
  save = FALSE
  filetype = "png"
  path = NULL
  xlab = NULL
  biomass_ylab = 'Biomass',
  cpue_ylab = 'CPUE'

}
