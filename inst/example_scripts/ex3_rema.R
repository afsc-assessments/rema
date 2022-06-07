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
