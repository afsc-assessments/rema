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

# m1: zeros as NAs ----
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

# m2: small constant = 0.0001, cv = 0.5 ----
input2 <- prepare_rema_input(model_name = 'small constant = 0.0001, cv = 0.5',
                            biomass_dat = nonsst,
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 0.5)))
m2 <- fit_rema(input2)
check_convergence(m2)
m2sum <- tidy_rema(m2)

# m3: small constant = 0.1, cv = 3 ----
input3 <- prepare_rema_input(model_name = 'small constant = 0.1, cv = 3',
                            biomass_dat = nonsst,
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.1, 3)))
m3 <- fit_rema(input3)
check_convergence(m3)
m3sum <- tidy_rema(m3)
compare <- compare_rema_models(rema_models = list(m1, m2, m3))
compare$plots$total_predicted_biomass

# m4: tweedie ----
input4 <- prepare_rema_input(model_name = 'tweedie',
                             biomass_dat = nonsst,
                             zeros = list(assumption = 'tweedie'))
input4$map
# input4 <- prepare_rema_input(model_name = 'tweedie',
#                              biomass_dat = nonsst,
#                              zeros = list(assumption = 'tweedie',
#                                           options_tweedie = list(zeros_cv = 1.2)))

# zeros = list(assumption = 'tweedie',
#              options_tweedie = list(fix_pars = c(1))))

m4 <- fit_rema(input4)
check_convergence(m4)
m4sum <- tidy_rema(m4)
m4sum$parameter_estimates
plot_rema(m4sum)$biomass_by_strata

compare <- compare_rema_models(rema_models = list(m1, m2, m3, m4))
compare$plots$total_predicted_biomass

# in log-space
bind_rows(m1sum$biomass_by_strata, m2sum$biomass_by_strata,
                     m3sum$biomass_by_strata, m4sum$biomass_by_strata) %>%
  dplyr::mutate(log_pred_lci = log_pred - 1.96 * sd_log_pred,
         log_pred_uci = log_pred + 1.96 * sd_log_pred) %>%
  ggplot(aes(x = year, y = log_pred, col = model_name, fill = model_name,  ymin = log_pred_lci, ymax = log_pred_uci)) +
  geom_ribbon(col = NA, alpha = 0.25) +
  geom_line() +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::scale_colour_viridis_d()

