# Exploring alternative assumptions for zero biomass observations.
# Example: other rockfish (all species except shortspine thornyhead) in the Bering Sea

# Example R scripts and data files:
# inst/example_scripts
# inst/example_data

# set up ----
library(rema)
library(dplyr)
library(cowplot) # provides helpful plotting utilities
library(ggplot2)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# create directory for analysis
# out_path <- "test_rema"
if(!exists("out_path")) out_path = getwd()
if(!dir.exists(out_path)) dir.create(out_path)

# copy all data files to working directory
rema_path <- find.package('rema')
example_data_files <- list.files(path = file.path(rema_path, "example_data"))
example_data_files
file.copy(from = file.path(path = file.path(rema_path, "example_data"),
                           example_data_files),
          to = file.path(file.path(out_path), example_data_files),
          overwrite = TRUE)

setwd(out_path)

# read data ----

biomass_dat <- read.csv(paste0("bsai_shortraker_biomass.csv"))
cpue_dat <- read.csv(paste0("bsai_shortraker_cpue.csv"))

# BTS only ----

input <- prepare_rema_input(model_name = 'BSAI shortraker BTS',
                            biomass_dat = biomass_dat,
                            zeros = list(assumption = 'NA'))
names(input)

# (2) Fit REMA model
m <- fit_rema(input)

# (3) Get tidied data.frames from the REMA model output
output <- tidy_rema(rema_model = m)
names(output)
output$parameter_estimates # estimated fixed effects parameters
output$biomass_by_strata # data.frame of predicted and observed biomass by stratum
output$total_predicted_biomass # total predicted biomass (same as biomass_by_strata for univariate models)

# (4) Generate model plots
plots <- plot_rema(tidy_rema = output,
                   # optional y-axis label
                   biomass_ylab = 'Biomass (t)')
plots$biomass_by_strata

# (5) Get one-step-ahead (OSA) residuals
osa <- get_osa_residuals(m)
osa$residuals
osa$residuals$biomass %>% filter(is.nan(residual))

cowplot::plot_grid(osa$plots$biomass_resids,
                   osa$plots$biomass_qqplot,
                   osa$plots$biomass_hist,
                   osa$plots$biomass_fitted)

osa <- get_osa_residuals(m, options = list(method = "oneStepGeneric"))
osa <- get_osa_residuals(m, options = list(method = "fullGaussian"))
osa <- get_osa_residuals(m, options = list(method = "oneStepGaussianOffMode"))
osa <- get_osa_residuals(m, options = list(method = "oneStepGaussian"))

# AI shortraker ----

admb_re <- read_admb_re(filename = 'aisr_rwout.rep',
                        # optional label for the single biomass survey stratum
                        biomass_strata_names = 'Aleutians Islands',
                        model_name = 'ADMB: AI shortraker')

input <- prepare_rema_input(model_name = 'M1: Base', admb_re = admb_re)
m <- fit_rema(input)
output <- tidy_rema(rema_model = m)
plots <- plot_rema(tidy_rema = output, biomass_ylab = 'Biomass (t)')
plots$biomass_by_strata

osa <- get_osa_residuals(m)
osa$residuals
osa$residuals$biomass %>% filter(is.nan(residual))

cowplot::plot_grid(osa$plots$biomass_resids,
                   osa$plots$biomass_qqplot,
                   osa$plots$biomass_hist,
                   osa$plots$biomass_fitted)

osa$residuals$biomass %>% filter(is.nan(residual))
osa <- get_osa_residuals(m, options = list(method = "oneStepGeneric"))
osa <- get_osa_residuals(m, options = list(method = "fullGaussian"))
osa <- get_osa_residuals(m, options = list(method = "oneStepGaussianOffMode"))
osa <- get_osa_residuals(m, options = list(method = "oneStepGaussian"))

