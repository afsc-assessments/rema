
# Test examples of running REMA using existing ADMB RE model output report files
# (rwout.rep)

# Example R scripts and data files:
# inst/example_scripts
# inst/example_data

library(rema)
# install.packages('cowplot')
library(cowplot) # provides helpful plotting utilities

ggplot2::theme_set(cowplot::theme_cowplot() + cowplot::background_grid())

# Ex 1 RE ----

# Univariate version of the random effects model (i.e., single survey, single
# stratum) using an existing ADMB rwout.rep file. Example: Aleutian Islands
# shortraker (aisr.rep) with NMFS bottom trawl survey estimates

# (1) Read in existing rwout.rep files, which is the report file generated from
# the ADMB version of the random effects model
?read_admb_re
admb_re <- read_admb_re(filename = 'inst/example_data/aisr_rwout.rep',
                      # optional label for the single biomass survey stratum
                      biomass_strata_names = 'Aleutians Islands')
names(admb_re)

# (2) Prepare REMA model inputs
?prepare_rema_input # note alternative methods for bringing in survey data observations
input <- prepare_rema_input(model_name = 'ai_shortraker_re',
                            admb_re = admb_re)
names(input)

# (3) Fit REMA model
?fit_rema
m <- fit_rema(input)
names(m)

# (4) Check convergence criteria if you so wish
?check_convergence
check_convergence(m)

# (5) Get tidied data.frames from the REMA model output
?tidy_rema
output <- tidy_rema(rema_model = m)
names(output)
output$parameter_estimates # estimated fixed effects parameters
output$biomass_by_strata # data.frame of predicted and observed biomass by stratum
output$total_predicted_biomass # total predicted biomass (same as biomass_by_strata for univariate models)

# (6) Generate model plots
?plot_rema
plots <- plot_rema(tidy_rema = output,
                   # optional y-axis label
                   biomass_ylab = 'Biomass (t)')

plots$biomass_by_strata
plots$total_predicted_biomass
plots$biomass_by_cpue_strata

# Ex 2 REM ----

# Multivariate version of the random effects model (REM) with a single survey
# and multiple strata. Example using Bering Sea and Aleutian Islands shortspine
# thornyhead

admb_re <- read_admb_re(filename = 'inst/example_data/bsaisst_rwout.rep',
                      biomass_strata_names = c('AI survey', 'EBS slope survey', 'S. Bering Sea (AI survey)'))

input <- prepare_rema_input(model_name = 'bsai_sst_rem',
                            admb_re = admb_re)

m <- fit_rema(input)
check_convergence(m)

output <- tidy_rema(rema_model = m)
output$parameter_estimates # estimated fixed effects parameters
output$biomass_by_strata # data.frame of predicted and observed biomass by stratum
output$total_predicted_biomass

plots <- plot_rema(tidy_rema = output,
                   # optional y-axis label
                   biomass_ylab = 'Biomass (t)')

plots$biomass_by_strata
plots$biomass_by_strata + ggplot2::facet_wrap(~strata, ncol = 1, scales = 'free_y')
plots$total_predicted_biomass
plots$total_predicted_biomass + ggplot2::ggtitle('BSAI Shortspine thornyhead predicted biomass')

# Ex 3 REMA ----

# Multi-survey and multi-strata version of the random effects model (REMA).
# Example using GOA shortraker rockfish, which uses the same strata definitions
# for the biomass and CPUE survey.
admb_re <- read_admb_re(filename = 'inst/example_data/goasr_rwout.rep',
                        biomass_strata_names = c('CGOA', 'EGOA', 'WGOA'),
                        cpue_strata_names = c('CGOA', 'EGOA', 'WGOA'))

input <- prepare_rema_input(model_name = 'GOA shortraker',
                            multi_survey = 1,
                            admb_re = admb_re,
                            sum_cpue_index = TRUE,
                            # one process error parameters (log_PE) estimated
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)),
                            # three scaling parameters (log_q) estimated, indexed as
                            # follows for each biomass survey stratum:
                            q_options = list(pointer_q_biomass = c(1, 2, 3)))

m <- fit_rema(input)
check_convergence(m)

output <- tidy_rema(m)
output$parameter_estimates

plots <- plot_rema(output, biomass_ylab = 'Biomass (t)', cpue_ylab = 'Relative Population Weights')
plots$biomass_by_strata
plots$cpue_by_strata
cowplot::plot_grid(plots$biomass_by_strata,
                   plots$cpue_by_strata,
                   ncol = 1)

plots$total_predicted_biomass
plots$total_predicted_cpue
plots$biomass_by_cpue_strata

# Ex 4 REMA ----

# Multi-survey and multi-strata version of the random effects model (REMA).
# Example using GOA shortspine thornyhead, which has different strata
# definitions for the biomass and CPUE surveys.
admb_re <- read_admb_re(filename = 'inst/example_data/goasst_rwout.rep',
                      biomass_strata_names = c('CGOA (0-500 m)', 'CGOA (501-700 m)', 'CGOA (701-1000 m)',
                                               'EGOA (0-500 m)', 'EGOA (501-700 m)', 'EGOA (701-1000 m)',
                                               'WGOA (0-500 m)', 'WGOA (501-700 m)', 'WGOA (701-1000 m)'),
                      cpue_strata_names = c('CGOA', 'EGOA', 'WGOA'))
admb_re$biomass_dat
length(unique(admb_re$biomass_dat$strata))
length(unique(admb_re$cpue_dat$strata))
input <- prepare_rema_input(model_name = 'GOA shortspine thornyhead',
                            multi_survey = 1,
                            admb_re = admb_re,
                            sum_cpue_index = TRUE,
                            # three process error parameters (log_PE) estimated, indexed
                            # as follows for each biomass survey stratum:
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            # three scaling parameters (log_q) estimated, indexed as
                            # follows for each biomass survey stratum:
                            q_options = list(pointer_q_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                                             pointer_q_cpue = c(1, 2, 3)))

m <- fit_rema(input)
check_convergence(m)
output <- tidy_rema(m)
output$parameter_estimates
plots <- plot_rema(output, biomass_ylab = 'Biomass (t)', cpue_ylab = 'Relative Population Weight')
plots$biomass_by_strata
plots$cpue_by_strata
plots$biomass_by_cpue_strata
plots$total_predicted_biomass
plots$total_predicted_cpue

cowplot::plot_grid(plots$biomass_by_strata + facet_wrap(~strata, nrow = 1),
                   plots$cpue_by_strata, nrow = 2)
cowplot::plot_grid(plots$biomass_by_cpue_strata, plots$cpue_by_strata, nrow = 2)
