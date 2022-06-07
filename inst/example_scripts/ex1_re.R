# Simple univariate version of the random effects model using an existing ADMB
# rwout.rep file. Example: Aleutian Islands shortraker (aisr.rep) with NMFS bottom
# trawl survey estimates

# inst/example_scripts
# inst/example_data

#
re_dat <- read_re_dat(filename = 'inst/example_data/aisr.rep')

prepare_rema_input(model_name = 'AI_shortraker_RE',
                   re_dat = re_dat)

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

