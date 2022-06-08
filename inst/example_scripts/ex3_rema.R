# REMA example
library(rema)
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
input$data$biomass_obs
input$data$wt_biomass
input$data$wt_cpue

# test values - remove when fxn is complete
q_options = list(pointer_q_cpue = c(1, 2, 3),
                 pointer_q_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                 initial_pars = NULL,
                 fix_pars = NULL,
                 penalty_options = 'normal_prior',
                 penalty_values = c(c(1.0, 0.08), c(1.5, 0.08), c(1.0, 0.8)))

# test values - remove when fxn is complete
PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                  penalty_options = 'normal_prior',
                  penalty_values = c(c(1.0, 0.08), c(1.5, 0.08), c(1.0, 0.8)))
