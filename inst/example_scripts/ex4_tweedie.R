# Tweedie for observation error

# This script explores the following:
# 1) Observation error confidence intervals for Tweedie
# 2) The effect of the zero biomass CV assumption (defaults to 1.5)
# 2) Other examples, including REs without zeros, REM, and REMAs (start by using
# SEO YE)

# Example R scripts and data files:
# inst/example_scripts
# inst/example_data

library(rema)
library(ggplot2)
# install.packages('cowplot')
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

# simple RE ----

nonsst <- read.csv('ebsshelf_orox.csv')

input <- prepare_rema_input(model_name = 'tweedie zeros_cv = 1.5',
                            biomass_dat = nonsst,
                            zeros =
                              list(assumption = 'tweedie'
                                   # use the following arguments to change the
                                   # assumed cv on zeros (zeros_cv), fix
                                   # logit_tweedie_p (fix_pars), and change the
                                   # starting value or fixed values of
                                   # logit_tweedie_p (initial_pars):

                                   # , options_tweedie =
                                   #   list(zeros_cv = 2.0,
                                   #        fix_pars = 1,
                                   #        # if tweedie p = 1.9, #
                                   #        log_tweedie_p = log((1.9-1)/(2-1.9))
                                   #        initial_pars = log((1.9 - 1) / (2 - 1.9)))
                              ))

m <- fit_rema(input)
check_convergence(m)
msum <- tidy_rema(m)
msum$parameter_estimates
plot_rema(msum)$biomass_by_strata

# assumption about tweedie CIs on
newci <- msum$biomass_by_strata %>%
  select(year, obs) %>%
  bind_cols(data.frame(tweedie_sd = m$report()$biomass_sd)) %>%
  mutate(tweedie_lci = obs - tweedie_sd * 1.96,
         tweedie_uci = obs + tweedie_sd * 1.96) %>%
  mutate(tweedie_lci = ifelse(tweedie_lci < 0, 0, tweedie_lci))

ci2 <- plot_rema(msum)$total_predicted_biomass +
  geom_point(data = newci, aes(x = year, y = obs)) +
  geom_errorbar(data = newci,
                aes(x = year, y = obs,
                    ymin = tweedie_lci,
                    ymax = tweedie_uci)) +
  ggtitle('Tweedie observation error')

ci1 <- plot_rema(msum)$biomass_by_strata +
  ggtitle('Lognormal observation error')

plot_grid(ci1, ci2, ncol = 1)

# effect of zeros_cv ----

input <- prepare_rema_input(model_name = 'zeros_cv = 0.5',
                            biomass_dat = nonsst,
                            zeros =  list(assumption = 'tweedie',
                                          options_tweedie = list(zeros_cv = 0.5)))
input$map$logit_tweedie_p; input$par$logit_tweedie_p
cbind(input$data$biomass_obs, input$data$biomass_cv)

m1 <- fit_rema(input)
check_convergence(m1)
msum1 <- tidy_rema(m1)
msum1$parameter_estimates
input$map$logit_tweedie_p

input <- prepare_rema_input(model_name = 'zeros_cv = 1.4',
                            biomass_dat = nonsst,
                            zeros =  list(assumption = 'tweedie',
                                          options_tweedie = list(zeros_cv = 1.4)))
m2 <- fit_rema(input)
check_convergence(m2)
msum2 <- tidy_rema(m2)

input <- prepare_rema_input(model_name = 'zeros_cv = 1.5',
                            biomass_dat = nonsst,
                            zeros =  list(assumption = 'tweedie',
                                          options_tweedie = list(zeros_cv = 1.5)))
m3 <- fit_rema(input)
check_convergence(m3)
msum3 <- tidy_rema(m3)

input <- prepare_rema_input(model_name = 'zeros_cv = 2.0',
                            biomass_dat = nonsst,
                            zeros =  list(assumption = 'tweedie',
                                          options_tweedie = list(zeros_cv = 2.0)))
m4 <- fit_rema(input)
check_convergence(m4)
msum4 <- tidy_rema(m4)

input <- prepare_rema_input(model_name = 'zeros_cv = 5.0',
                            biomass_dat = nonsst,
                            zeros =  list(assumption = 'tweedie',
                                          options_tweedie = list(zeros_cv = 5.0)))
m5 <- fit_rema(input)
check_convergence(m5)
msum5 <- tidy_rema(m5)
msum5$parameter_estimates

compare <- compare_rema_models(rema_models = list(m2, m3, m4, m5))

# ex2 YE ----

ye_biomass <- read.csv('seo_ye_biomass.csv')

input <- prepare_rema_input(model_name = 'YE tweedie',
                            biomass_dat = ye_biomass %>%
                              filter(strata %in% c('CSEO', 'EYKT')) %>%
                              mutate(cv = 0.5),
                            zeros =  list(assumption = 'tweedie'))#,
                                          # options_tweedie = list(fix_pars = c(1))))
input$map$logit_tweedie_p; input$par$logit_tweedie_p
cbind(input$data$biomass_obs, input$data$biomass_cv)

m1 <- fit_rema(input)
check_convergence(m1)
m1$runtime
msum1 <- tidy_rema(m1)
msum1$parameter_estimates
input$map$logit_tweedie_p
plot_rema(msum1)$biomass_by_strata

input <- prepare_rema_input(model_name = 'YE lognorm',
                            biomass_dat = ye_biomass)
input$map$logit_tweedie_p; input$par$logit_tweedie_p
cbind(input$data$biomass_obs, input$data$biomass_cv)

m2 <- fit_rema(input)
check_convergence(m2)
msum2 <- tidy_rema(m2)
msum2$parameter_estimates

compare <- compare_rema_models(rema_models = list(m1, m2))
compare$plots$biomass_by_strata
