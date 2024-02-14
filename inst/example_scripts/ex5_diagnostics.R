# currently on the diagnostics branch!

# resources and references:
# https://github.com/NOAA-FIMS/TMB_training/tree/main
# https://github.com/eco4cast/Statistical-Methods-Seminar-Series/tree/main/auger-methe-ssm-in-tmb

library(readr)
library(dplyr)
library(ggplot2)
library(TMB)
# library(rema)
library(cowplot)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# dyn.unload(dynlib(here::here('src', 'rema')))
# TMB::compile(here::here('src', 'rema.cpp'))
# dyn.load(dynlib(here::here('src', 'rema')))
biomass_dat <- read_csv('inst/example_data/goa_sst_biomass.csv')
cpue_dat <- read_csv('inst/example_data/goa_sst_rpw.csv')

# library(rema)
unique(biomass_dat$strata)
input <- prepare_rema_input(model_name = 'GOA thornyhead',
                            biomass_dat = biomass_dat,
                            # shared process error SD across all strata
                            # PE_options = list(pointer_PE_biomass = c(rep(1,9)))
                            # shared process error SD within regions
                            PE_options = list(pointer_PE_biomass = c(rep(1,3), rep(2,3), rep(3,3)))
                            )

m1 <- fit_rema(input)
out1 <- tidy_rema(m1)
param <- out1$parameter_estimates %>% dplyr::mutate(region = c('CGOA', 'EGOA', 'WGOA')) # always lists the strata in alphabetical order!
param
plots1 <- plot_rema(tidy_rema = out1)
plots1$biomass_by_strata
plots1$total_predicted_biomass

# Simulation self-check ----

# 1) fit model to real data
# 2) use parameters fitted model to simulate data (TMB is awesome because it
# does that for us! well, almost)
# 3) use that simulated data to re-estimate parameters
# 4) check for bias between the two sets of parameters estimates (can we recover
# our parameters?)
set.seed(415)

# this is going to be a headache to generalize different strata, multi survey
# version, tweedie, extra CV. we can discuss.
for(i in 1:200) {
  sim <- m1$simulate(complete = TRUE)
  names(sim)
  newinput <- input
  tmp <- matrix(data = exp(sim$log_biomass_obs), ncol = ncol(input$data$biomass_obs), nrow = nrow(input$data$biomass_obs))
  colnames(tmp) <- colnames(input$data$biomass_obs)
  newinput$data$biomass_obs <- tmp
  # is there a reason we'd want to output and plot the simulated data?
  newm1 <- fit_rema(newinput)
  newout1 <- tidy_rema(newm1)
  if(i == 1) {
    simout <- newout1$parameter_estimates %>%
      dplyr::mutate(region = c('CGOA', 'EGOA', 'WGOA'),
                    sim = i)
  } else {
    simout <- simout %>%
      dplyr::bind_rows(newout1$parameter_estimates %>%
                         dplyr::mutate(region = c('CGOA', 'EGOA', 'WGOA'),
                                       sim = i))
  }
}

ggplot(simout, aes(x = estimate)) +
  geom_histogram(fill = 'white', col = 'black') +
  geom_vline(data = param, aes(xintercept = estimate), col = 'red', size = 2) +
  facet_wrap(~region, scales = 'free_x')

# Notes on simulations:
#
# > I'm not too sure about the posterior predictive checks or other uses for
# simulations in this particular case but if you want to dig into them more that
# would be great
#
# > How many simulations do we need to test bias? I don't know. Also can we get
# rid of all the dplyr junk printed to the screen?
#
# > Output bias statistic e.g. relative error?

# Identifiability and Uncertainty -----

# Parameters are "identifiable" if they have a solution (e.g., a minimum) and
# aren't highly correlated with another parameter

# A note about TMB::tmbprofile: Because parameters have the same name, need to
# use this "lincomb" argument. The line below is a likelihood profile (and
# confidence interval, etc.) for the second stratum's log_PE
profile <- TMB::tmbprofile(m1, lincomb = c(0,1,0))
plot(profile)
confint(profile)
# ^ note this is an alternative way of calculating variance to the standard
# method (see below) which relies on the delta method and assumes asymptotic
# approximation. the likelihood profile does not make the asymptotic normal
# assumption. the routine refits the model at a series of fixed values for the
# parameter of interest. the confidence interval is determined by the quantile
# of the chi-square distribution with 1 degree of free at p=0.95 (divided by 2
# because 2*log-likelihood difference asymptotically has this distribution). in
# this way it is also a good diagnostic to test whether the asymptotic normal
# assumption is valid or not (in this case it is)
m1$sdrep
# sdrep will give you the parameters at the scale (e.g., log or logit) at which
# they are estimated), whereas tidy_rema(m1)$parameter_estimates will give the
# user the estimates transformed to the arithmetic scale
tidy_rema(m1)$parameter_estimates

# Checking the Laplace approximation ----

# TMB allows the user to test whether the Laplace approximation of the marginal
# log-likelihood and joint log-likelihood is ok. TMB::checkConsistency simulates
# data and calculates the gradients for each of the simulated data sets. If the
# average gradient is approximately zero, then the Laplace approximation is
# assumed to be ok. A chi-square test for gradient bias is performed for the
# marginal and joint log-likelihood and estimates of bias are provided for
# fixed effects parameters. Increasing sample sizes will increase power of the
# chi-square test for bias.
check <- TMB::checkConsistency(m1)
check
# Its very concerning to me that we are unable to invert the information matrix
# (same as the Hessian here right?!)
summary(check)

# One-step-ahead residuals ----

# I've already built a wrapper function (not to say it couldn't use some work!) to do this get_osa_residuals() but
# consistently have issues with residual patterns... I guess they may in fact
# reflect true model misspecification!
resids <- get_osa_residuals(m1)  # uses cdf methods
resids$residuals
# check out the NaNs! in this case its always the first observation. this isn't
# always the case but usually.
resids$residuals$biomass %>% filter(is.nan(residual))
cowplot::plot_grid(resids$plots$biomass_resids, resids$plots$biomass_qqplot, ncol = 1)

# One thought i've had is that these issues have to due with our initial
# conditions. To test this we can try fixing the first biomass observation using
# the map object in TMB to see if it appropriately deals with our OSA problems
# and consistency in the Laplace approximation
newinput <- input
newinput$model_name <- "GOA thornyhead (fixed init)"
nstrata <- ncol(newinput$par$log_biomass_pred) # index for the number of biomass strata
nyr <- nrow(newinput$par$log_biomass_pred) # index for the number of years in the model
tmp <- as.factor(c(rep(NA, nstrata), # map off (aka fix) the first year of biomass
                   1:(length(newinput$map$log_biomass_pred) - nstrata))) # numbers in the map are uniquely estimated parameters
tmp <- matrix(data = tmp, ncol = ncol(newinput$par$log_biomass_pred), nrow = nyr,
              byrow = TRUE) # must be byrow! otherwise it won't be indexed properly
newinput$map$log_biomass_pred <- as.factor(as.vector(tmp))
newm1 <- fit_rema(newinput)
newout1 <- tidy_rema(newm1)
newout1$parameter_estimates
param

newresids <- get_osa_residuals(newm1)  # uses cdf methods
newresids$residuals
# sure enough, at least in this case it got rid of the NaNs.... yikes!
newresids$residuals$biomass %>% filter(is.nan(residual))
cowplot::plot_grid(newresids$plots$biomass_resids, newresids$plots$biomass_qqplot, ncol = 1)

TMB::checkConsistency(newm1) # but this still fails, so that's good :)

compare <- compare_rema_models(list(m1, newm1))
compare$plots$total_predicted_biomass # gives identical results from a management perspective
compare$plots$biomass_by_strata + facet_wrap(~strata, scales = 'free_y')
