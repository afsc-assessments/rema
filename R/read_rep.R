# Read ADMB .rep file and return an R object of type 'list'
#
# Code from Steve Martell, D'Arcy N. Webbe

# called by read_re_dat.R
# fn = name of ADMB output file to be read
# returns object of type "list" with ADMB outputs therein

read_rep <- function(fn) {
  options(warn = -1) # Suppress the NA message in the coercion to double
  repfile <- scan(fn, what = "character", flush = TRUE, blank.lines.skip = FALSE, quiet = TRUE, na.strings = c("nan","-nan"))
  inan <- which(is.na(repfile)) # Identify any nan entries so that they are not picked up as objects
  idx <- sapply(as.double(repfile), is.na)
  idx[inan] <- FALSE
  vnam <- repfile[idx] # list names
  nv <- length(vnam) # number of objects
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
