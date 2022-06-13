
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
library(ggplot2)
# install.packages('cowplot')
library(cowplot) # provides helpful plotting utilities like plot_grid() and nice ggplot2 themes

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# create directory for analysis: e.g., out_path <- "/path/to/save/output"
if(!exists("out_path")) out_path = getwd()
if(!dir.exists(out_path)) dir.create(out_path)
setwd(out_path)

# copy all data files to working directory
rema_path <- find.package("rema")

example_data_files <- list.files(path = file.path(rema_path, "example_data"))
file.copy(from = file.path(path = file.path(rema_path, "example_data"),
                           example_data_files),
          to = file.path(out_path),
          overwrite = TRUE)

# confirm you are in the working directory and it has the the example rwout.rep
# and csv files
list.files()

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
# each stratum - 8 fixed effects parameters
input1 <- prepare_rema_input(model_name = 'REMA_SEO_YE_strataPE',
                            multi_survey = TRUE,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat)
m1 <- fit_rema(input1)
check_convergence(m1)
output1 <- tidy_rema(m1)
output1$parameter_estimates
plots1 <- plot_rema(output1, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
plots1$biomass_by_strata
plots1$cpue_by_strata
cowplot::plot_grid(plots1$biomass_by_strata + facet_wrap(~strata, ncol = 1),
                   plots1$cpue_by_strata + facet_wrap(~strata, ncol = 1),
                   nrow = 1)
plots1$total_predicted_cpue # note that total cpue is not available because nominal cpue is not summable

# Model 2: one process error shared across all strata and separate CPUE survey scaling parameters for
# each stratum - 5 fixed effects parameters
input2 <- prepare_rema_input(model_name = 'REMA_SEO_YE_singlePE',
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
plots2 <- plot_rema(output2, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
cowplot::plot_grid(plots2$biomass_by_strata + facet_wrap(~strata, nrow = 1),
                   plots2$cpue_by_strata + facet_wrap(~strata, nrow = 1),
                   nrow = 2)

# Model 3: Separate process errors for northern and southern distributs (i.e.
# one shared process error in EYKT and NSEO and one shared process error in CSEO
# and NSEO) - 6 fixed effects parameters
colnames(input2$data$biomass_obs) # must index strata in alphabetical order
input3 <- prepare_rema_input(model_name = 'REMA_SEO_YE_northsouthPE',
                             multi_survey = TRUE,
                             biomass_dat = biomass_dat,
                             cpue_dat = cpue_dat,
                             # pointer_PE_biomass is an index that assigns a
                             # process error parameter to each biomass survey
                             # stratum (MUST INDEX STRATA NAMES BY ALPHABETICAL
                             # ORDER! "CSEO" "EYKT" "NSEO" "SSEO") double check
                             # order using colnames(input2$data$biomass_obs)
                             PE_options = list(pointer_PE_biomass = c(1, 2, 2, 1)))
m3 <- fit_rema(input3)
check_convergence(m3)
output3 <- tidy_rema(m3)
output3$parameter_estimates
plots3 <- plot_rema(output3, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
cowplot::plot_grid(plots3$biomass_by_strata + facet_wrap(~strata, nrow = 1),
                   plots3$cpue_by_strata + facet_wrap(~strata, nrow = 1),
                   nrow = 2)

# (3) Compare models
compare <- compare_rema_models(rema_models = list(m1, m2, m3),
                               biomass_ylab = 'ROV biomass',
                               cpue_ylab = 'IPHC setline survey CPUE')

compare$aic # The single process error model has the lowest AIC and is the simplest model
names(compare$plots)

compare$plots$total_predicted_biomass

cowplot::plot_grid(compare$plots$biomass_by_strata +
                     facet_wrap(~strata, nrow = 1) +
                     theme(legend.position = 'top'),
                   compare$plots$cpue_by_strata +
                     facet_wrap(~strata, nrow = 1) +
                     theme(legend.position = 'none'),
                   nrow = 2) # use , rel_heights = c(0.52, 0.48)) if you want to get anal-retentive about it...

# (4) Fit  models that don't include the IPHC survey CPUE

# This model will have 1 PE per strata - 4 fixed effects parameters
input4 <- prepare_rema_input(model_name = 'REMA_SEO_YE_noCPUE_strataPE',
                             multi_survey = FALSE, # <--- removes CPUE
                             biomass_dat = biomass_dat)
m4 <- fit_rema(input4)
check_convergence(m4)
output4 <- tidy_rema(m4)
output4$parameter_estimates
plots4 <- plot_rema(output4, biomass_ylab = 'ROV biomass')
plots4$biomass_by_strata

# This model will have one shared PE - 1 fixed effects parameter
input5 <- prepare_rema_input(model_name = 'REMA_SEO_YE_noCPUE_singlePE',
                             multi_survey = FALSE, # <--- removes CPUE
                             biomass_dat = biomass_dat,
                             PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)) # <-- shared process error
                             )
m5 <- fit_rema(input5)
check_convergence(m5)
output5 <- tidy_rema(m5)
output5$parameter_estimates
plots5 <- plot_rema(output5, biomass_ylab = 'ROV biomass')
plots5$biomass_by_strata

# (5) model comparison

compare <- compare_rema_models(rema_models = list(m1, m2, m3, m4, m5))
compare$aic # no AIC because models were fit to different data
compare$plots$biomass_by_strata
compare$plots$total_predicted_biomass
compare$plots$cpue_by_strata # does not exist...check output:
compare$output$cpue_by_strata
compare_rema_models(rema_models = list(m4, m5))$aic

