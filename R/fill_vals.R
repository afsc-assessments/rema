# function to fill list component with a factor. used for TMB map generation
fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}
