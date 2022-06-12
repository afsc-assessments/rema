
# Test example of reading biomass and survey data from .csv, fitting multiple
# models, and doing model comparison with AIC. Example uses yelloweye rockfish
# data in Southeast Alaska Outside (SEO) waters. Multi-survey, multi-strata
# (REMA) fit to ROV survey biomass estimates and IPHC setline survey CPUE
# (numbers per hook), with four strata (EYKT = East Yakutat, NSEO = Northern
# SEO, CSEO = Central SEO, SSE = Southern SEO)

# Example R scripts and data files:
# inst/example_scripts
# inst/example_data

library(rema)
# install.packages('cowplot')
library(cowplot) # provides helpful plotting utilities like plot_grid() and nice ggplot2 themes

ggplot2::theme_set(cowplot::theme_cowplot() + cowplot::background_grid())

# (1) Read the biomass and cpue survey data from file. See ?prepare_rema_input
# for information on required columns
biomass_dat <- read.csv('inst/example_data/seo_ye_biomass.csv')
str(biomass_dat)
unique(biomass_dat$strata)
cpue_dat <- read.csv('inst/example_data/seo_ye_cpue.csv')
str(cpue_dat)
unique(cpue_dat$strata)

# (2) Fit alternative models

# Model 1: separate process errors and CPUE survey scaling parameters for
# each stratum (8 fixed effects parameters)
input1 <- prepare_rema_input(model_name = 'SEO_Yelloweye',
                            multi_survey = TRUE,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat)
m1 <- fit_rema(input1)
check_convergence(m1)
output1 <- tidy_rema(m1)
output1$parameter_estimates
plots1 <- plot_rema(output1, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey cpue')
plots1$biomass_by_strata
plots1$cpue_by_strata
cowplot::plot_grid(plots1$biomass_by_strata + facet_wrap(~strata, nrow = 1),
                   plots1$cpue_by_strata + facet_wrap(~strata, nrow = 1),
                   nrow = 2)
plots1$total_predicted_cpue # note that total cpue is not available because nominal cpue is not summable

# Model 2: one process error shared across all strata and separate CPUE survey scaling parameters for
# each stratum (5 fixed effects parameters)
input2 <- prepare_rema_input(model_name = 'SEO_Yelloweye',
                            multi_survey = TRUE,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            # pointer_PE_biomass is an index that assigns a
                            # process error parameter to each biomass survey
                            # stratum (default = one PE for each stratum, as in
                            # Model 1)
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)))
m2 <- fit_rema(input2)
check_convergence(m2)
output2 <- tidy_rema(m2)
output2$parameter_estimates
plots2 <- plot_rema(output2, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey cpue')
cowplot::plot_grid(plots2$biomass_by_strata + facet_wrap(~strata, nrow = 1),
                   plots2$cpue_by_strata + facet_wrap(~strata, nrow = 1),
                   nrow = 2)

# model comparison
