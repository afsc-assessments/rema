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
library(tidyverse)
library(devtools)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# dyn.unload(dynlib(here::here('src', 'rema')))
TMB::compile(here::here('src', 'rema.cpp'))
dyn.load(dynlib(here::here('src', 'rema')))
biomass_dat <- read_csv('inst/example_data/goa_sst_biomass.csv')
# biomass_dat <- biomass_dat %>% filter(strata == "CGOA (0-500 m)")
cpue_dat <- read_csv('inst/example_data/goa_sst_rpw.csv')

# library(rema)
unique(biomass_dat$strata)
devtools::document()
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
m1$sdrep
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

# made quiet for temporary case at least
# make data frame to hold things in
re_est <- matrix(NA, 200, 3); colnames(re_est) <- c("CGOA", "EGOA", "WGOA")

# do the things
suppressMessages(for(i in 1:200) {
  sim <- m1$simulate(complete = TRUE)
  # names(sim)
  newinput <- input
  tmp <- matrix(data = exp(sim$log_biomass_obs), ncol = ncol(input$data$biomass_obs), nrow = nrow(input$data$biomass_obs))
  colnames(tmp) <- colnames(input$data$biomass_obs)
  newinput$data$biomass_obs <- tmp
  # is there a reason we'd want to output and plot the simulated data?
  newm1 <- fit_rema(newinput, do.sdrep = F) # need to run with something to supress report? do.sdrep = F fine for running but then freaks out below
  re_est[i, ] <- newm1$env$last.par[1:3]
  # newout1 <- tidy_rema(newm1)
  # if(i == 1) {
  #   simout <- newout1$parameter_estimates %>%
  #     dplyr::mutate(region = c('CGOA', 'EGOA', 'WGOA'),
  #                   sim = i)
  # } else {
  #   simout <- simout %>%
  #     dplyr::bind_rows(newout1$parameter_estimates %>%
  #                        dplyr::mutate(region = c('CGOA', 'EGOA', 'WGOA'),
  #                                      sim = i))
  # }
})

# ggplot everything
re_est_long <- as.data.frame(re_est) %>% pivot_longer(cols = 1:3, names_to = "region", values_to = "est_PE")
ggplot(re_est_long, aes(exp(est_PE))) + # simout, aes(x = estimate)) +
  geom_histogram(fill = 'white', col = 'black') +
  geom_vline(data = param, aes(xintercept = estimate), col = 'red', size = 2) +
  facet_wrap(~region, scales = 'free_x')

# Notes on simulations:
#
# > I'm not too sure about the posterior predictive checks or other uses for
# simulations in this particular case but if you want to dig into them more that
# would be great

# so looking at marie's code for posterior predictives, it seems pretty similar?
# posterior predictive is basically checking that the simulated observations are working
# asking: if we have those parameters, do we get the output (observations) we expect?
# so one step less than above tbh?
# note two options: one where par_est is fixed and one where par_est is drawn from distribution

# prep
nrep <- 5000 # simulate time series nrep times
# make matrix
mat_mean <- matrix(NA, nrep, 9) # 9 strata: means will go here
mat_sd <- matrix(NA, nrep, 9) # and sds will go here
# get real means & sds of logged observations (i.e., convert to model log world and then take summary stats)
true_vals <- biomass_dat %>% group_by(strata) %>%
  summarise(true_mean = mean(log(biomass), na.rm = TRUE), true_sd = sd(log(biomass), na.rm = TRUE))

# loop
for (i in 1:nrep) {

  # draw parameters: if drawing (like bayesian posterior predictive)
  par_est <- m1$env$last.par # set up with factor structure and all other parameters
  par_est[1:3] <- log(rnorm(3, mean = out1$parameter_estimates$estimate, sd = out1$parameter_estimates$std_err)) # get parameter estimates into par_est: note that they are logged within model and exponential in the estimates

  # so we want to basically specify the estimates we'll use to simulate
  # two options... not sure which is best here (see notes above)
  # if drawing
  sim <- m1$simulate(par=par_est) # this tells m1 to simulate the data using the estimated model parameters (log_PE estimates)
  # if not drawing:
  # sim <- m1$simulate(complete = TRUE) # this uses the last.par values: so uses "real" states??
  # note: "The default parameter values used for the simulation is obj$env$last.par" https://kaskr.github.io/adcomp/Simulation.html
  # not quite sure that using the par_est is fair in this case, since it preserves the starting values? not sure it matters either? (doesn't really make a huge difference with either simulation option tbh)

  # summary stats into the matrix
  # note simulation output is in LOG world: this keeps things normally distributed?
  # put means in the mean matrix
  mat_mean[i, ] <- apply(sim$log_biomass_obs, 2, mean, na.rm = TRUE)
  # put sd in the sd matrix
  mat_sd[i, ] <- apply(sim$log_biomass_obs, 2, sd, na.rm = TRUE)

  # could be cool to do a bunch of light lines and then the real distribution for output plot?

}

# make things into a data frame/ggplot situation
df_mean <- as.data.frame(mat_mean); colnames(df_mean) <- unique(biomass_dat$strata)
df_sd <- as.data.frame(mat_sd); colnames(df_sd) <- unique(biomass_dat$strata)
# rotate things... rip
df_mean_long <- df_mean %>% pivot_longer(cols = 1:9, names_to = "strata", values_to = "mean_log_biomass")
df_sd_long <- df_sd %>% pivot_longer(cols = 1:9, names_to = "strata", values_to = "sd_log_biomass")

# plots
ggplot(df_mean_long, aes(mean_log_biomass)) + geom_histogram(fill = 'white', col = 'black') +
  facet_wrap(~strata) +
  geom_vline(data = true_vals, aes(xintercept = true_mean), col = "red", size = 2) +
  theme_classic() + labs(x = "mean log biomass over time period")
ggplot(df_sd_long, aes(sd_log_biomass)) + geom_histogram(fill = 'white', col = 'black') +
  facet_wrap(~strata) +
  geom_vline(data = true_vals, aes(xintercept = true_sd), col = "red", size = 2) +
  theme_classic() + labs(x = "sd log biomass over time period") # its walking too much???

# relative error calc?
mean_relerror <- apply((mat_mean - matrix(rep(t(true_vals$true_mean)), nrow = nrep, ncol = 9, byrow = T))/matrix(rep(t(true_vals$true_mean)), nrow = nrep, ncol = 9, byrow = T), 2, mean)
sd_relerror <- apply((mat_sd - matrix(rep(t(true_vals$true_sd)), nrow = nrep, ncol = 9, byrow = T))/matrix(rep(t(true_vals$true_sd)), nrow = nrep, ncol = 9, byrow = T), 2, mean)

# > How many simulations do we need to test bias? I don't know. Also can we get
# rid of all the dplyr junk printed to the screen?
# --> not sure how many simulations. posterior check runs really fast w/5000 simuations (<10 seconds)
# --> used suppressMessages to remove dplyr things, suppressed TMB output info in a funky way but hopefully not horrible way?
#
# > Output bias statistic e.g. relative error?
# --> right now just the plots plus simple relative error calculation for each strata/output

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
# note hessian issue:
m1$hessian <- TRUE
opt <- do.call("optim", m1)
opt
opt$hessian ## <-- FD hessian from optim
m1$he()    ## <-- Analytical hessian

# using MCMC approach: tmbstan :)
# note need stan chain set up, RIP
library(tmbstan)
# okay so can put fitted model in? "All methods provided by the rstan package can be applied to a fitted object"
# https://cran.r-project.org/web/packages/tmbstan/tmbstan.pdf
fit.la <- tmbstan(obj = m1, chains = 3, init = m1$par, laplace = T)
fit.mcmc <- tmbstan(obj = m1, chains = 3, init = m1$par, laplace = F) # doesn't fit super great but post warm up looks okay?
# can check chains, trace plot
rstan::traceplot(fit.la, pars=names(m1$par), inc_warmup=TRUE) # should be squiggly catipliers
rstan::traceplot(fit.mcmc, pars=names(m1$par), inc_warmup=TRUE) # should be squiggly catipliers
# get the parameter draws
post.la <- as.data.frame(fit.la)[, c(1:3)]; post.la$sim_type <- "LA"
post.mcmc <- as.data.frame(fit.mcmc)[, c(1:3)]; post.mcmc$sim_type <- "MCMC"
plot(post.la$`log_PE[1]`, post.la$lp__) # side note: check identifiability here
# qq comparisions
# basically, want a to compare the quantiles from each of the parameters between the two fits
qv <- seq(from = 0, to = 1, by = 0.01)
post.quantiles <- data.frame(parameter = c(rep(1, length(qv)), rep(2, length(qv)), rep(3, length(qv))),
                             quantile.vals = c(rep(qv, 3)),
                             mcmc = c(quantile(post.mcmc$`log_PE[1]`, probs = qv),
                                      quantile(post.mcmc$`log_PE[2]`, probs = qv),
                                      quantile(post.mcmc$`log_PE[3]`, probs = qv)),
                             lpax = c(quantile(post.la$`log_PE[1]`, probs = qv),
                                      quantile(post.la$`log_PE[2]`, probs = qv),
                                      quantile(post.la$`log_PE[3]`, probs = qv)))
ggplot(post.quantiles, aes(mcmc, lpax, col = as.factor(parameter))) +
  facet_wrap(~parameter, scales = "free") + geom_abline(intercept = 0, slope = 1, col = "lightgray") + geom_point() +
  theme_classic() + theme(legend.position = "none") + scale_color_manual(values = c("lightpink", "lightblue", "lightgreen"))
# so looks good, but this isn't really a metric, it's a plot?
# second plot will be the comparision of the pairs plots
post.la.long <- post.la %>% pivot_longer(!sim_type, names_to = "parameter", values_to = "value") %>% arrange(parameter)
post.mcmc.long <- post.mcmc %>% pivot_longer(!sim_type, names_to = "parameter", values_to = "value") %>% arrange(parameter)
# 3000 samples each --> rotate some thing around
post.la.long$parameter_comp <- c(post.la.long$parameter[3001:9000], post.la.long$parameter[1:3000])
post.la.long$value_comp <- c(post.la.long$value[3001:9000], post.la.long$value[1:3000])
post.mcmc.long$parameter_comp <- c(post.mcmc.long$parameter[3001:9000], post.mcmc.long$parameter[1:3000])
post.mcmc.long$value_comp <- c(post.mcmc.long$value[3001:9000], post.mcmc.long$value[1:3000])
post.all.long <- rbind(post.mcmc.long, post.la.long) # combine everyone
# get the means
post.la.means <- post.la.long %>% group_by(sim_type, parameter, parameter_comp) %>% summarise(value = mean(value), value_comp = mean(value_comp))
post.mcmc.means <- post.mcmc.long %>% group_by(sim_type, parameter, parameter_comp) %>% summarise(value = mean(value), value_comp = mean(value_comp))
post.means <- rbind(post.la.means, post.mcmc.means)
# make the plot
ggplot(post.all.long, aes(value, value_comp, col = sim_type)) +
  geom_point(pch = 1, size = 0.5) + facet_grid(rows = vars(as.factor(parameter)), cols = vars(as.factor(parameter_comp))) +
  geom_point(data = post.means, size = 1.2, pch = 1, color = "black") +
  theme_classic() + scale_color_manual(values = c(alpha("lightpink", 0.5), alpha("lightblue", 0.5)))
# or use GGally dependency for plot
library(GGally)
GGally::ggpairs(rbind(post.la, post.mcmc), aes(color = sim_type)) +
  scale_color_manual(values = c(alpha("lightpink", 0.8), alpha("lightblue", 0.8))) +
  scale_fill_manual(values = c(alpha("lightpink", 0.8), alpha("lightblue", 0.8)))
# so again, this is a plot that kind of works (emphasis on kind of), but not a metric....
bias_stat <- rbind(apply(post.mcmc[1:3], 2, mean), apply(post.la[1:3], 2, mean))
bias_stat_tidy <- (bias_stat[1, ] - bias_stat[2, ]) # need to clean up
# should also be able to get a p-value for if the normal distributions are different from each other... t test?
t.test(post.la$`log_PE[1]`, post.mcmc$`log_PE[1]`)
t.test(post.la$`log_PE[2]`, post.mcmc$`log_PE[2]`)
t.test(post.la$`log_PE[3]`, post.mcmc$`log_PE[3]`)
# is there a way to summarise across all three? and is this the right stat? bc really just saying the means are the same (which is also what the means/bias table says)

# One-step-ahead residuals ----

# Methods for calculating OSA residuals in REMA have been implemented in TMB
# through the `TMB::oneStepPredict()` function with the `fullGaussian` method
# function. The use of OSA residuals, sometimes referred to as forecast
# residuals or prediction errors, is crucial for validation of state-space
# models like REMA (Thygesen et al. 2017). Instead of comparing the observed and
# expected values at the same time step, OSA residuals use forecasted values
# based on all previous observations (i.e. they exclude the observed value at
# the current time step from prediction). Traditional residuals (e.g. Pearson's
# residuals) are inappropriate for state-space models, because process error
# variance may be over-estimated in cases where the model is mispecified, thus
# leading to artificially small residuals (see Section 3 of Thygesen et al. 2017
# for an example).

# A wrapper function is available in rema to get leverages TMB::oneStepPredict()
# using the fullGaussian method as default
resids <- get_osa_residuals(m1)
resids$residuals$biomass # residuals! :)

# A quantile-quantile (Q-Q) plot can help determine if residuals are likely to
# have come from a normal distribution, which is a key assumption in the model.
resids$plots$qq

# An autocorrelation function (ACF) plot of residuals can show if there is
# autocorrelation in the residuals of time series data when fitting a regression
# model. Autocorrelation can occur when the value of a variable in the current
# time period is similar to its value in previous periods. If the ACF plot of
# residuals shows significant lags, it could mean that the estimated model
# violates the assumption of no autocorrelation in the errors.
resids$plots$acf
resids$plots$histo
resids$plots$biomass_resids
resids$plots$biomass_fitted

# One could also compare other residual methods
resids2 <- get_osa_residuals(m1, options = list("cdf"))
resids2$residuals$biomass$residual
# when the normality assumptions underlying the fullGaussian methods are met,
# they should yield identical results to cdf methods
plot(resids$residuals$biomass$residual, resids2$residuals$biomass$residual)
cor(x=resids$residuals$biomass$residual[!is.na(resids$residuals$biomass$residual)], y=resids2$residuals$biomass$residual[!is.na(resids2$residuals$biomass$residual)])

# OSA simulation test ----

# to do with simple & complex case: (1) simulate 1000 cases, refit them & check
# validations, e.g., residuals (2) make matrix of residuals: row replicate, col
# data point (3) KS test for normality, distribution of KS test p values should
# be uniform enough
r <- TMB::oneStepPredict(obj = m1, observation.name = "obsvec",
                         data.term.indicator = "keep",
                         discrete = FALSE, method = "fullGaussian")
nresid <- length(r$residual)
nsim = 1000
simout <- matrix(NA, nrow = nsim, ncol = nresid)

# suppressMessages(
for(i in 1:nsim) {
  sim <- m1$simulate(complete = TRUE)
  # sim$biomass_obs[1,]
  # exp(sim$log_biomass_obs)[1,]
  newinput <- input
  tmp <- matrix(data = exp(sim$log_biomass_obs), ncol = ncol(input$data$biomass_obs), nrow = nrow(input$data$biomass_obs))
  colnames(tmp) <- colnames(input$data$biomass_obs)
  newinput$data$biomass_obs <- tmp
  newinput$data$obsvec <- sim$log_biomass_obs[!is.na(sim$log_biomass_obs)]
  newm1 <- fit_rema(newinput, do.sdrep = FALSE)
  tmpr <- TMB::oneStepPredict(obj = newm1, observation.name = "obsvec",
                           data.term.indicator = "keep",
                           discrete = FALSE, method = "fullGaussian")
  simout[i,] <- tmpr$residual
}
# )

ksout <- vector(length = nsim)
# OSA residuals should be distributed iid N(0,1) under a model with valid
# assumptions

# one-sided KS test for normality
# H0: the distribution of the residuals in normal
# Ha: the distribution of the residuals in not normal
# If p > 0.05 there is no evidence to reject the null hypothesis that the
# distribution of the data is normal.
for(i in 1:nsim) {
  kstmp <- ks.test(x=simout[i,], "pnorm", mean=0, sd=1)
  ksout[i] <- kstmp$p.value
}

ks.test(ksout, "runif", 0, 1)
length(ksout[ksout>=0.05])/length(ksout)

saveRDS(list(simout = simout, ksout = ksout),
        file="inst/example_scripts/ex5_osa_sim.RData")
