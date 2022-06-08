# Simple univariate version of the random effects model using an existing ADMB
# rwout.rep file. Example: Aleutian Islands shortraker (aisr.rep) with NMFS bottom
# trawl survey estimates

# inst/example_scripts
# inst/example_data

#
library(rema)

# univariate ----
re_dat <- read_re_dat(filename = 'inst/example_data/aisr.rep')

input <- prepare_rema_input(model_name = 'AI_shortraker_RE',
                   re_dat = re_dat)

m <- fit_rema(input)

# tidy_rema <- function(m$input)
m$sdrep$value

names(m)
m$par
m$is_sdrep
m$sdrep
obj <- TMB::MakeADFun(data = input$data,
                      parameters = input$par,
                      random = input$random,
                      dll = 'rema',
                      map = input$map)
opt <- with(obj, nlminb(par, fn, gr))
adrep <- TMB::sdreport(obj)
adrep$value
obj$report()

# REM ----

re_dat <- read_re_dat(filename = 'inst/example_data/bsaisst.rep')

input <- prepare_rema_input(model_name = 'BSAI_shortspine_thornyhead',
                            re_dat = re_dat)

obj <- TMB::MakeADFun(data = input$data,
                      parameters = input$par,
                      random = input$random,
                      dll = 'rema',
                      map = input$map)
opt <- with(obj, nlminb(par, fn, gr))
adrep <- TMB::sdreport(obj)
adrep
adrep$sd
adrep$value
adrep$par.random # biomass random effects

best <- obj$env$last.par.best # maximum likelihood estimates
adrep$gradient.fixed
names(adrep)
adrep$value[names(adrep$value) == 'biomass_pred']
adrep$value

# REMA ----

re_dat <- read_re_dat(filename = 'inst/example_data/goasst.rep')
re_dat$biomass_dat
length(unique(re_dat$biomass_dat$strata))
length(unique(re_dat$cpue_dat$strata))
input <- prepare_rema_input(model_name = 'GOA shortspine thornyhead',
                            multi_survey = 1,
                            # three process error parameters (log_PE) estimated, indexed
                            # as follows for each biomass survey stratum:
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            # three scaling parameters (log_q) estimated, indexed as
                            # follows for each biomass survey stratum:
                            q_options = list(pointer_q_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            re_dat = re_dat)

obj <- TMB::MakeADFun(data = input$data,
                      parameters = input$par,
                      random = input$random,
                      dll = 'rema',
                      map = input$map)
opt <- with(obj, nlminb(par, fn, gr))
adrep <- TMB::sdreport(obj)
adrep
adrep$sd
adrep$value
adrep$par.random # biomass random effects
