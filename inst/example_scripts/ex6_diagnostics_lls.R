# Diagnostic tools for the REMA model
#
# The purpose of this script is to explore and develop model
# diagnostic functions for rema users

# **currently on the diagnostics branch**

# resources and references:
# https://github.com/NOAA-FIMS/TMB_training/tree/main
# https://github.com/eco4cast/Statistical-Methods-Seminar-Series/tree/main/auger-methe-ssm-in-tmb

library(readr)
library(dplyr)
library(ggplot2)
library(TMB)
library(rema)
library(cowplot)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# compile model locally
# dyn.unload(dynlib(here::here('src', 'rema')))
TMB::compile(here::here('src', 'rema.cpp'))
dyn.load(dynlib(here::here('src', 'rema')))

biomass_dat <- read_csv('inst/example_data/goa_sst_biomass.csv')
cpue_dat <- read_csv('inst/example_data/goa_sst_rpw.csv')

unique(biomass_dat$strata)
input <- prepare_rema_input(model_name = 'GOA thornyhead',
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            multi_survey = TRUE,
                            # Process error variation alternation configurations:
                            # (1) shared process error SD across all strata
                            # PE_options = list(pointer_PE_biomass = c(rep(1,9)))
                            # (2) shared process error SD within an area
                            PE_options = list(pointer_PE_biomass = c(rep(1,3), rep(2,3), rep(3,3))),
                            # (3) one process error SD for 0-500 shared across
                            # areas, and one process error SD for > 500 m shared
                            # across all areas (the idea being that younger at
                            # shallower depths may be more variable than older
                            # fish at deeper depths)
                            # PE_options = list(pointer_PE_biomass = rep(c(1,2,2),3)),
                            # LLS scaling parameter: reminder that that the bts
                            # and lls are not stratified in the same way... so
                            # this pointer allows the user to map biomass strata
                            # to cpue strata
                            q_options = list(pointer_biomass_cpue_strata = c(rep(1,3), rep(2,3), rep(3,3))),
                            # estimate additional obs error on the LLS
                            extra_cpue_cv = list(assumption = 'extra_cv'),
                            # rpns and rpws are summable
                            sum_cpue_index = TRUE
                            )

m1 <- fit_rema(input)
out1 <- tidy_rema(m1)

# note on this "region_input": not sure if this could be gleaned from the
# pointers in the prepare_rema_input or an argument to tidy_rema... TBD
region_input <- c(rep(c('CGOA', 'EGOA', 'WGOA'), 2), "COMBINED") # always lists the strata in alphabetical order!
param <- out1$parameter_estimates %>% dplyr::mutate(region = region_input)
param
plots1 <- plot_rema(tidy_rema = out1)
plots1$biomass_by_strata + facet_wrap(~strata, scales = 'free_y', dir = 'v')
plots1$cpue_by_strata + facet_wrap(~strata, scales = 'free_y')

# plot the additional estimated observation error on the LLS
out1cv <- tidy_extra_cv(out1)
plot_extra_cv(out1cv)$cpue_by_strata

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
for(i in 1:150) {
  sim <- m1$simulate(complete = TRUE)
  names(sim)
  newinput <- input
  tmp <- matrix(data = exp(sim$log_biomass_obs), ncol = ncol(input$data$biomass_obs), nrow = nrow(input$data$biomass_obs))
  colnames(tmp) <- colnames(input$data$biomass_obs)
  newinput$data$biomass_obs <- tmp
  tmp <- matrix(data = exp(sim$log_cpue_obs), ncol = ncol(input$data$cpue_obs), nrow = nrow(input$data$cpue_obs))
  colnames(tmp) <- colnames(input$data$cpue_obs)
  newinput$data$cpue_obs <- tmp
  # is there a reason we'd want to output and plot the simulated data?
  newm1 <- fit_rema(newinput)
  newout1 <- tidy_rema(newm1)
  if(i == 1) {
    simout <- newout1$parameter_estimates %>% dplyr::mutate(region = region_input, sim = i)
  } else {
    simout <- simout %>% dplyr::bind_rows(newout1$parameter_estimates %>% dplyr::mutate(region = region_input, sim = i))
  }
}

simout <- simout %>% dplyr::mutate(name = paste0(region, "_", parameter))
param <- param %>% dplyr::mutate(name = paste0(region, "_", parameter))

ggplot(simout, aes(x = estimate)) +
  geom_histogram(fill = 'white', col = 'black') +
  geom_vline(data = param, aes(xintercept = estimate), col = 'red', size = 2) +
  facet_wrap(~name, scales = 'free_x')

# Notes on simulations:
#
# > Pursue posterior predictive checks. unclear if the self-check is a good
# diagnostic tool (maybe just better for checking code for errors or poor
# convergence?)... there did not seem to be good agreement within the ssmesa
# group about the utility of the simulation self-check as a model validation
# tool
#
# > Add flags for simulating process error vs. observations in the cpp file
#
# > How many simulations do we need to test bias?
#
# > Get rid of all the dplyr junk printed to the screen
#
# > Output bias statistic e.g., relative error?

# Identifiability and Uncertainty -----

# Parameters are "identifiable" if they have a solution (e.g., a minimum) and
# aren't highly correlated with another parameter

# A note about TMB::tmbprofile: Because parameters have the same name, need to
# use this "lincomb" argument. The line below is a likelihood profile (and
# confidence interval, etc.) for the second stratum's log_PE
profile <- TMB::tmbprofile(m1, lincomb = c(0,1,0,0,0,0,0))
plot(profile);confint(profile)
profile <- TMB::tmbprofile(m1, lincomb = c(0,0,0,0,0,0,1))
plot(profile);confint(profile)
# # can also look at the difference between parameters - Anders N says don't do this!
# profile <- TMB::tmbprofile(m1, lincomb = c(0,0,1,0,0,-1,0))
plot(profile);confint(profile)
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

# Flagging this paper, which outlines another potential method we can use to
# test the accuracy of the Laplace approximation:
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197954#pone-0197954-g001

# TMB allows the user to test accuracy of the Laplace approximation of the marginal
# log-likelihood and joint log-likelihood. TMB::checkConsistency simulates
# data and calculates the gradients for each of the simulated data sets. If the
# average gradient is approximately zero, then the Laplace approximation is
# assumed to be ok. A chi-square test for gradient bias is performed for the
# marginal and joint log-likelihood and estimates of bias are provided for
# fixed effects parameters. Increasing sample sizes will increase power of the
# chi-square test for bias.
check <- TMB::checkConsistency(m1)
check
# In this case the simulation does not appear to be correct!! We tried turning
# off SIMULATE block on random effects and that caused the entire thing to fail.

# in other case I'm not able to invert the information matrix
# (same as the Hessian here right?)
summary(check)

# One-step-ahead residuals ----

TMB::oneStepPredict(obj = m1,
                    observation.name = "obsvec",
                    data.term.indicator = "keep",
                    # discrete = FALSE,
                    method = "cdf")


# TO DO:
# > Change default in OSA function to 'generic' based on Kasper's recommendation
# (via Cole)
#
# I've built a wrapper function (not to say it couldn't use some work!) to do this get_osa_residuals() but
# consistently have issues with residual patterns... I guess they may in fact
# reflect true model misspecification!
resids <- get_osa_residuals(m1, options = list(method = "cdf"))  # uses cdf methods
resids$residuals
# check out the NaNs! in this case its always the first observation. this isn't
# always the case but usually.
resids$residuals$biomass %>% filter(is.nan(residual))
cowplot::plot_grid(resids$plots$biomass_resids, resids$plots$biomass_qqplot, ncol = 1)
cowplot::plot_grid(resids$plots$cpue_resids, resids$plots$cpue_qqplot, ncol = 1)

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
cowplot::plot_grid(resids$plots$cpue_resids, resids$plots$cpue_qqplot, ncol = 1)

TMB::checkConsistency(newm1) # but this still fails, so that's good :)

compare <- compare_rema_models(list(m1, newm1))
compare$plots$total_predicted_biomass # gives identical results from a management perspective
compare$plots$biomass_by_strata + facet_wrap(~strata, scales = 'free_y')
compare$plots$cpue_by_strata + facet_wrap(~strata, scales = 'free_y')

