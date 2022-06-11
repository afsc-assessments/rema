# Simple univariate version of the random effects model using an existing ADMB
# rwout.rep file. Example: Aleutian Islands shortraker (aisr.rep) with NMFS bottom
# trawl survey estimates

# inst/example_scripts
# inst/example_data

#
library(rema)

re_dat <- read_re_dat(filename = 'inst/example_data/aisr.rep')

input <- prepare_rema_input(model_name = 'AI_shortraker_RE',
                            re_dat = re_dat)

m <- fit_rema(input)
check_convergence(m)

m_output <- tidy_rema(rema_model = m)

m_plots <- plot_rema(tidy_output = m_output,
                     biomass_ylab = 'Biomass (t)')
plot_rema <- function()
  # plotting function
  biomass_ylab = 'Biomass (t)'
  cpue_ylab = 'Relative Population Weight'

  # how many plot rows should I have?
  p1 <- ggplot2::ggplot(data = biomass_summary,
                        ggplot2::aes(x = year, y = pred)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pred_lci, ymax = pred_uci),
                         col = 'grey', fill = 'grey') +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~strata, nrow = NULL) +
    ggplot2::geom_point(ggplot2::aes(x = year, y = obs)) +
    ggplot2::geom_errorbar(ggplot2::aes(x = year, ymin = obs_lci, ymax = obs_uci)) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(x = NULL, title = biomass_index_description) +
    theme_minimal_grid()
  p1


  p2 <- ggplot2::ggplot(data = cpue_summary,
                        ggplot2::aes(x = year, y = estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci),
                         col = 'grey', fill = 'grey') +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~strata, nrow = 1) +
    ggplot2::geom_point(ggplot2::aes(x = year, y = cpue)) +
    ggplot2::geom_errorbar(ggplot2::aes(x = year, ymin = obslci, ymax = obsuci)) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(x = NULL, y = 'CPUE index')
  p2


  library(cowplot)

  cowplot::plot_grid(p1, p2, nrow = 2)

m$sdrep$value
names(m)
m$rep
names(m)
m$par
m$is_sdrep
m$sdrep
obj <- TMB::MakeADFun(data = input$data,
                      parameters = input$par,
                      random = input$random,
                      dll = 'rema',
                      map = input$map)
opt <- with(obj, nlminb(par, fn, gr))
adrep <- TMB::sdreport(obj)
adrep$value
obj$report()

# BSAI shortspine thornyhead REM ----

re_dat <- read_re_dat(filename = 'inst/example_data/bsaisst.rep')

input <- prepare_rema_input(model_name = 'BSAI_shortspine_thornyhead',
                            re_dat = re_dat)


# GOA shortspine thornyhead REMA ----
re_dat <- read_re_dat(filename = 'inst/example_data/goasst.rep',
                      biomass_strata_names = c('CGOA1', 'CGOA2', 'CGOA3',
                                               'EGOA1', 'EGOA2', 'EGOA3',
                                               'WGOA1', 'WGOA2', 'WGOA3'),
                      cpue_strata_names = c('CGOA', 'EGOA', 'WGOA'))
re_dat$biomass_dat
length(unique(re_dat$biomass_dat$strata))
length(unique(re_dat$cpue_dat$strata))
input <- prepare_rema_input(model_name = 'GOA shortspine thornyhead',
                            multi_survey = 1,
                            re_dat = re_dat,
                            sum_cpue_index = TRUE,
                            # three process error parameters (log_PE) estimated, indexed
                            # as follows for each biomass survey stratum:
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            # three scaling parameters (log_q) estimated, indexed as
                            # follows for each biomass survey stratum:
                            q_options = list(pointer_q_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3))
)

m3 <- fit_rema(input)
check_convergence(m3)

# REMA GOA shortraker ----
re_dat <- read_re_dat(filename = 'inst/example_data/goasr.rep',
                      biomass_strata_names = c('CGOA', 'EGOA', 'WGOA'),
                      cpue_strata_names = c('CGOA', 'EGOA', 'WGOA'))
re_dat$biomass_dat
length(unique(re_dat$biomass_dat$strata))
length(unique(re_dat$cpue_dat$strata))
input <- prepare_rema_input(model_name = 'GOA shortraker',
                            multi_survey = 1,
                            re_dat = re_dat,
                            sum_cpue_index = TRUE,
                            # one process error parameters (log_PE) estimated
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)),
                            # three scaling parameters (log_q) estimated, indexed as
                            # follows for each biomass survey stratum:
                            q_options = list(pointer_q_biomass = c(1, 2, 3))
)

m4 <- fit_rema(input)
check_convergence(m4)

