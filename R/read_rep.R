#' Read ADMB .rep file and return an R object of type 'list'
#'
#' Code modified from original function provided by Steve Martell, D'Arcy N. Webber

#' called by \code{\link{read_re_dat}}
#'
#' @param fn full path and name of ADMB output file to be read
#'
#' @return object of type "list" with ADMB outputs therein
#' @export
#'
#' @examples
#' \dontrun{
#' read_rep(fn = 'inst/example_data/goasr.rep')
#' }
read_rep <- function(fn) {
  options(warn = -1) # Suppress the NA message in the coercion to double
  repfile <- scan(fn, what = "character", flush = TRUE, blank.lines.skip = FALSE, quiet = TRUE, na.strings = c("nan","-nan"))
  inan <- which(is.na(repfile)) # Identify any nan entries so that they are not picked up as objects
  idx <- sapply(as.double(repfile), is.na)
  idx[inan] <- FALSE
  vnam <- repfile[idx] # list names
  nv <- length(vnam) # number of objects

  # errors specific to rema
  if(nv != length(unique(vnam))) {
    dups <- vnam[duplicated(vnam)]
    stop(paste0("The following variable(s) are duplicated in the rwout.rep file provided by the user: ", dups, ". The duplicate entries must be removed in order for read_re_dat() to function properly. Please check the re.tpl for a duplicate write_R() statement to fix future rwout.rep files."))
  }

  A <- list()
  ir <- 0
  for (i in 1:nv) {
    ir <- match(vnam[i], repfile)
    if (i != nv) {
      irr <- match(vnam[i+1], repfile)
    } else {
      irr <- length(repfile) + 1 # next row
    }
    dum <- NA
    # vector/scalar objects
    if (irr-ir == 2) {
      dum <- as.double(scan(fn, skip = ir, nlines = 1, quiet = TRUE, what = ""))
    }
    # matrix objects
    if (irr-ir > 2) {
      dum <- as.matrix(read.table(fn, skip = ir, nrow = irr-ir-1, fill = TRUE, row.names = NULL))
    }
    if (is.numeric(dum)) { # Logical test to ensure dealing with numbers
      A[[vnam[i]]] <- dum
    }
  }
  options(warn = 0)
  return(A)
}
