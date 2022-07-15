# Exploring alternative assumptions for zero biomass observations.
# Example: other rockfish (all species except shortspine thornyhead) in the Bering Sea

# Example R scripts and data files:
# inst/example_scripts
# inst/example_data

# set up ----
library(rema)
library(dplyr)
library(cowplot) # provides helpful plotting utilities

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
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

# read in data ----
nonsst <- read.csv('ebsshelf_orox.csv')

# fit REMA models ----

# model 1: zeros as NAs
input1 <- prepare_rema_input(model_name = 'zeros as NAs',
                             biomass_dat = nonsst,
                             zeros = list(assumption = 'NA'))
m1 <- fit_rema(input1)
check_convergence(m1)
tidy_rema(m1)
m1$report()
m1$sdrep
m1sum <- tidy_rema(m1)
plot_rema(m1sum)

# model 2: small constant = 0.0001, cv = 0.5
input2 <- prepare_rema_input(model_name = 'small constant = 0.0001, cv = 0.5',
                            biomass_dat = nonsst,
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 0.5)))
m2 <- fit_rema(input2)

# model 3: small constant = 0.1, cv = 3
input3 <- prepare_rema_input(model_name = 'small constant = 0.1, cv = 3',
                            biomass_dat = nonsst,
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.1, 3)))
m3 <- fit_rema(input3)

compare <- compare_rema_models(rema_models = list(m1, m2, m3))
compare$plots$total_predicted_biomass

# model 4: tweedie
input4 <- prepare_rema_input(model_name = 'tweedie',
                             biomass_dat = nonsst %>%
                               mutate(cv = ifelse(biomass == 0, 1.5, cv)),
                             zeros = list(assumption = 'tweedie'))

input4 <- prepare_rema_input(model_name = 'tweedie',
                             biomass_dat = nonsst %>%
                               mutate(cv = ifelse(biomass == 0, 1.5, cv)),
                             zeros = list(assumption = 'tweedie',
                                          options_tweedie = list(fix_pars = c(1))))
m4 <- fit_rema(input4, do.fit = T)
m4 <- fit_rema(input4, do.fit = F)
names(m4)
m4$gr()
m4$fn()
m4$report()
cbind(input4$data$biomass_obs, m4$report()$biomass_sd)

