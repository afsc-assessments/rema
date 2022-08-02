# BSAI shortraker biomass estimation

# set up ----

# assessment year
YEAR <- 2022

libs <- c('readr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# install.packages("devtools")
# devtools::install_github("JaneSullivan-NOAA/rema", dependencies = TRUE)
library(rema)

# folder set up
dat_path <- paste0("inst/example_data/"); dir.create(dat_path)
# out_path <- paste0("results/", YEAR); dir.create(out_path)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# read data ----

biomass_dat <- read_csv(paste0(dat_path, "/bsai_shortraker_biomass.csv"))
cpue_dat <- read_csv(paste0(dat_path, "/bsai_shortraker_cpue.csv"))

# BTS only ----

# (1) Prepare REMA model inputs
?prepare_rema_input # note alternative methods for bringing in survey data observations
input <- prepare_rema_input(model_name = 'BSAI shortraker BTS',
                            biomass_dat = biomass_dat,
                            zeros = list(assumption = 'NA'))
names(input)

# (2) Fit REMA model
?fit_rema
m <- fit_rema(input)

# (3) Check convergence criteria if you so wish
?check_convergence
check_convergence(m)

# (4) Get tidied data.frames from the REMA model output
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

# (7) Get one-step-ahead (OSA) residuals
osa <- get_osa_residuals(m)
osa$residuals
cowplot::plot_grid(osa$plots$biomass_resids,
                   osa$plots$biomass_qqplot,
                   osa$plots$biomass_hist,
                   osa$plots$biomass_fitted)

# BTS + LLS ----

input <- prepare_rema_input(model_name = 'BSAI shortraker BTS + LLS',
                            multi_survey = 1, # fit to CPUE data? yes = 1)
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            sum_cpue_index = 1, # is the CPUE index summable (yes = 1, RPWs are summable)
                            zeros = list(assumption = 'NA'),
                            # sort(unique(biomass_dat$strata)) = "Central AI"
                            # "Eastern AI" "EBS Slope" "SBS" "Western AI". EBS
                            # Slope is stratum 3 for the biomass so we put it in
                            # the third position of the
                            # pointer_biomass_cpue_strata object and use NAs for
                            # the other 4 strata. It is the first (and only)
                            # stratum for the Longline survey, so we use the
                            # value of 1. See Details in ?prepare_rema_input
                            q_options = list(pointer_biomass_cpue_strata = c(NA, NA, 1, NA, NA)))

m2 <- fit_rema(input)
m2$report()
check_convergence(m2)

output <- tidy_rema(rema_model = m2)
output$parameter_estimates # estimated fixed effects parameters
output$biomass_by_strata # data.frame of predicted and observed biomass by stratum
output$total_predicted_biomass # total predicted biomass (same as biomass_by_strata for univariate models)

# (6) Generate model plots
?plot_rema
plots <- plot_rema(tidy_rema = output,
                   # optional y-axis label
                   biomass_ylab = 'Biomass (t)',
                   cpue_ylab = 'Relative population weights')
plots$biomass_by_strata
plots$cpue_by_strata

# (7) Get one-step-ahead (OSA) residuals
osa <- get_osa_residuals(m2)
osa$residuals
cowplot::plot_grid(osa$plots$biomass_resids,
                   osa$plots$biomass_qqplot,
                   osa$plots$biomass_hist,
                   osa$plots$biomass_fitted)
cowplot::plot_grid(osa$plots$cpue_resids,
                   osa$plots$cpue_qqplot,
                   osa$plots$cpue_hist,
                   osa$plots$cpue_fitted)

# Model comparisons -----

compare <- compare_rema_models(list(m,m2))
compare$plots$biomass_by_strata
compare$output$parameter_estimates
