% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_rema.R
\name{plot_rema}
\alias{plot_rema}
\title{Plot survey data and model output}
\usage{
plot_rema(
  tidy_rema,
  save = FALSE,
  filetype = "png",
  path = NULL,
  xlab = NULL,
  biomass_ylab = "Biomass",
  cpue_ylab = "CPUE"
)
}
\arguments{
\item{tidy_rema}{list out output from \code{\link{tidy_rema}}, which includes
model results but also inputs}

\item{save}{(optional) logical (T/F) save figures as \code{filetype} in
\code{path}. Default = FALSE. NOT YET IMPLEMENTED.}

\item{filetype}{(optional) character string; type of figure file. Default =
'png'. NOT YET IMPLEMENTED.}

\item{path}{(optional) directory path to location where figure files are to
be saved if \code{save = TRUE}. NOT YET IMPLEMENTED.}

\item{xlab}{(optional) label for x-axis of biomass and CPUE plots (e.g.
'Year'). Default = NULL.}

\item{biomass_ylab}{(optional) label for y-axis of biomass plots (e.g.
'Biomass (t)'). Default = 'Biomass'.}

\item{cpue_ylab}{(optional) label for y-axis of CPUE plots (e.g. 'Relative
Population Number'). Default = 'CPUE'.}
}
\value{
a list of ggplot2 plots or character string messages about the data.
  Except for parameter estimates, the objects output from
  \code{\link{tidy_rema}} are the same outputted from this function.
}
\description{
Takes list output from \code{\link{tidy_rema}} and returns a list of
\code{ggplot2} objects to be plotted or saved.
}
\examples{
\dontrun{
# placeholder for example
}
}
\seealso{
\code{\link{tidy_rema}}
}
