library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

# Packages for MCMC and diagnostics
library(tmbstan)
# tmbstan relies on rstan, which now needs to be installed through the R
# universe. if you are experiencing errors with MCMC, try uninstalling and then
# re-installing rstan using this code chunk (also don"t forget to restart your
# computer so that your path resets)
# install.packages("rstan", repos = c("https://stan-dev.r-universe.dev",
#                                     "https://cloud.r-project.org"))
library(rstan)
library(bayesplot)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 11) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

library(rema)

pcod_bio_dat <- read_csv("inst/example_data/ai_pcod_2022_biomass_dat.csv")
pcod_input <- prepare_rema_input(model_name = "p_cod",
                                 biomass_dat = pcod_bio_dat,
                                 # one strata
                                 PE_options = list(pointer_PE_biomass = 1)
)
# run the model
pcod_mod <- fit_rema(pcod_input)

# next, get the thornyhead data set up
thrn_bio_dat <- read_csv("inst/example_data/goa_thornyhead_2022_biomass_dat.csv")
thrn_cpue_dat <- read_csv("inst/example_data/goa_thornyhead_2022_cpue_dat.csv")
thrn_input <- prepare_rema_input(model_name = "thrnhead_rockfish",
                                 multi_survey = TRUE,
                                 biomass_dat = thrn_bio_dat,
                                 cpue_dat = thrn_cpue_dat,
                                 # RPWs are a summable/area-weighted effort index
                                 sum_cpue_index = TRUE,
                                 # three process error parameters (log_PE) estimated,
                                 # indexed as follows for each biomass survey stratum
                                 # (shared within an area across depths):
                                 PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                                 # scaling parameter options:
                                 q_options = list(
                                   # longline survey strata (n=3) indexed as follows for the
                                   # biomass strata (n=9)
                                   pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                                   # one scaling parameters (log_q) estimated, shared
                                   # over all three LLS strata
                                   pointer_q_cpue = c(1, 1, 1)),
                                 # estimate extra trawl survey observation error
                                 extra_biomass_cv = list(assumption = "extra_cv"),
                                 # estimate extra longline survey observation error
                                 extra_cpue_cv = list(assumption = "extra_cv")
)

# run the model
thrn_mod <- fit_rema(thrn_input)

# tidy output and plot fitted data
tidy_pcod <- tidy_rema(pcod_mod)
p1 <- plot_rema(tidy_pcod)$biomass_by_strata +
  ggtitle(label = "Model Fits to the AI Pcod Data",
          subtitle = "Trawl Survey Biomass Strata") +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
              col = 'goldenrod', fill = 'goldenrod', alpha = 0.4) +
  geom_line() +
  geom_point(aes(x = year, y = obs)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci))
p1
ggsave(paste0("vignettes/ex4_pcod_fits.png"), width = 6, height = 4, units = "in", bg = "white")

tidy_thrn <- tidy_rema(thrn_mod)
p2 <- plot_rema(tidy_thrn, biomass_ylab = 'Biomass (t)')$biomass_by_strata +
  ggtitle(label = "Model Fits to the GOA Thornyhead Data",
          subtitle = "Trawl Survey Biomass Strata") +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
              col = "#21918c", fill = "#21918c", alpha = 0.4) +
  geom_line() +
  geom_point(aes(x = year, y = obs)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci))
p2
p3 <- plot_rema(tidy_thrn, cpue_ylab = "RPW")$cpue_by_strata +
  ggtitle(label = NULL, subtitle = "Longline Survey Relative Population Weight (RPW) Strata") +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
              col = "#21918c", fill = "#21918c", alpha = 0.4) +
  geom_line() +
  geom_point(aes(x = year, y = obs)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci))
p3
cowplot::plot_grid(p2, p3, rel_heights = c(0.7, 0.3), ncol = 1)
ggsave(paste0("vignettes/ex4_thrn_fits.png"), width = 11, height = 10, units = "in", bg = "white")

# simulation -----

sim_test <- function(mod_name, replicates, cpue) {

  # storage things
  re_est <- matrix(NA, replicates, length(mod_name$par)) # parameter estimates
  grad <- matrix(NA, replicates, length(mod_name$par)) # gradients

  # go through the model
  suppressMessages(for(i in 1:replicates) {

    # i = 10; mod_name = thrn_mod; cpue = TRUE
    sim <- mod_name$simulate(complete = TRUE) # simulates the data

    # simulated biomass observations:
    tmp_biomass <- matrix(data = exp(sim$log_biomass_obs), ncol = ncol(mod_name$input$data$biomass_obs))
    colnames(tmp_biomass) <- colnames(mod_name$data$biomass_obs)
    # simulated cpue observations, when applicable:
    if (cpue) {tmp_cpue <- matrix(data = sim$cpue_obs, ncol = ncol(mod_name$input$data$cpue_obs))
    colnames(tmp_cpue) <- colnames(mod_name$data$cpue_obs)}
    # set up new data for input
    newinput <- mod_name$input
    newinput$data$biomass_obs <- tmp_biomass # biomass data
    if (cpue) {newinput$data$cpue_obs <- tmp_cpue} # cpue data

    # create "obsvec" which is used internally in cpp file as the observation
    # vector for all observation (log biomass + cpue) in the likelihood
    # functions (and required for OSA residuals). note the transpose t() needed
    # to get these matrices in the correct order (by row instead of by col) --
    # note obsvec is masked from users normally in prepare_rema_input()
    newinput$data$obsvec <- t(sim$log_biomass_obs)[!is.na(t(sim$log_biomass_obs))]
    if (cpue) {newinput$data$obsvec <- c(t(sim$log_biomass_obs)[!is.na(t(sim$log_biomass_obs))], t(sim$log_cpue_obs)[!is.na(t(sim$log_cpue_obs))])}

    # refit model
    mod_new <- fit_rema(newinput, do.sdrep = TRUE)

    # add parameter estimates to matrix
    if(mod_new$opt$convergence == 0) {
      re_est[i, ] <- mod_new$env$last.par.best[1:length(mod_name$par)]
    } else {
      re_est[i, ] <- rep(NA, length(mod_name$par))
    }
    grad[i, ] <- mod_new$final_gradient

  })

  re_est <- as.data.frame(re_est); re_est$type <- rep("recovered")
  grad <- as.data.frame(grad)

  return(list(est = re_est, grad = grad))

}

# run for pcod and prep data frame
# run simulation testing
pcod_sim <- sim_test(mod_name = pcod_mod, replicates = 500, cpue = FALSE)
par_ests <- pcod_sim$est; pcod_grad <- pcod_sim$grad
# get data frame with simulation values for each parameter
n_not_converged_pcod <- length(which(is.na(par_ests[,1]))); n_pcod <- length(par_ests[,1])
prop_converged_pcod <- 1-n_not_converged_pcod/n_pcod

mod_par_ests <- data.frame("log_PE1" = pcod_mod$env$last.par.best[1],
                           type = "model")
names(par_ests) <- names(mod_par_ests) # rename for ease
pcod_par_ests <- rbind(mod_par_ests, par_ests) # recovered and model in one data frame
pcod_par_ests$sp <- rep("AI Pcod")
names(pcod_grad) <- names(pcod_par_ests)[1]

# run for thorny and prep data frame -- same process
# run simulation testing
thrn_sim <- sim_test(mod_name = thrn_mod, replicates = 500, cpue = TRUE) # note some warnings
par_ests <- thrn_sim$est; thrn_grad <- thrn_sim$grad
# sometimes spits out: In stats::nlminb(model$par, model$fn, model$gr, control =
# list(iter.max = 1000,...: NA/NaN function evaluation
n_not_converged_thrn <- length(which(is.na(par_ests[,1]))); n_thrn <- length(par_ests[,1])
prop_converged_thrn <- 1-n_not_converged_thrn/n_thrn
nrow(par_ests)
# get data frame with simulation values for each parameter
mod_par_ests <- data.frame("log_PE1" = thrn_mod$env$last.par.best[1],
                           "log_PE2" = thrn_mod$env$last.par.best[2],
                           "log_PE3" = thrn_mod$env$last.par.best[3],
                           "log_q" = thrn_mod$env$last.par.best[4],
                           "log_tau_biomass" = thrn_mod$env$last.par.best[5],
                           "log_tau_cpue" = thrn_mod$env$last.par.best[6],
                           type = "model")
names(par_ests) <- names(mod_par_ests) # rename for ease
thrn_par_ests <- rbind(mod_par_ests, par_ests) # recovered and model parameters in one data frame
thrn_par_ests$sp <- rep("GOA Thornyhead")
names(thrn_grad) <- names(thrn_par_ests)[1:6]

# reogranize data
pcod_sim <- pcod_par_ests %>% pivot_longer(1, names_to = "parameter")
thrn_sim <- thrn_par_ests %>% pivot_longer(1:6, names_to = "parameter")
sim_dat <- rbind(pcod_sim, thrn_sim)
pcod_grad <- pcod_grad %>% pivot_longer(1, names_to = "parameter", values_to = "gradient") %>%
  mutate(high = ifelse(gradient > 1e-5, TRUE, FALSE))
thrn_grad <- thrn_grad %>% pivot_longer(1:6, names_to = "parameter", values_to = "gradient") %>%
  mutate(high = ifelse(gradient > 1e-35, TRUE, FALSE))
grad_dat <- rbind(pcod_grad,# %>% mutate(sp = "AI Pcod"),
                  thrn_grad) %>% # %>% mutate(sp = "GOA Thornyhead")) %>%
  select(-high,-parameter)# %>% mutate(type = 'recovered')
sim_dat %>%
  filter(type == "recovered") %>%
  bind_cols(grad_dat) %>%
  bind_rows(sim_dat %>% filter(type == "model")) %>%
  filter(is.na(value) | value < -9) %>%
  mutate(absgradient = abs(gradient)) %>%
  arrange(gradient) #%>% View


# Relative Error: ((om-em)/om)
sim_re <- sim_dat %>%
  pivot_wider(id_cols = c("sp", "parameter"),
              names_from = type, values_from = value) %>%
  unnest(cols = c(model, recovered)) %>%
  mutate(RE = (model-recovered)/model*100) %>%
  group_by(sp, parameter) %>%
  mutate(label = paste0("Median RE=", formatC(median(RE, na.rm = TRUE), format = "f", digits = 1), "%")) %>%
  suppressWarnings()

# plot distribution of simulated parameter estimates
plot_sim <- function(sim_dat, plot_title, fill_col = "#21918c") {
  ggplot(NULL, aes(parameter, value)) +
    # add distribution of recovered parameters
    geom_violin(data = sim_dat %>% filter(type == "recovered"),
                fill = fill_col, alpha = 0.6, draw_quantiles = 0.5) +
    # add "true values" from original model
    geom_point(data = sim_dat %>% filter(type == "model"),
               size = 2, col = "black") +
    facet_wrap(~ parameter, scales = "free", nrow = 1) +
    labs(x = NULL, y = "Parameter estimate", title = plot_title,
         subtitle = "Distribution of parameters estimates (median=horizontal line, true value=point)") +
    scale_x_discrete(labels = NULL, breaks = NULL)
}
# plot relative error
plot_re <- function(sim_re, fill_col = "#21918c") {
  ggplot(NULL, aes(parameter, RE)) +
    # add distribution of recovered parameters
    geom_boxplot(data = sim_re, alpha = 0.6, fill = fill_col,
                 na.rm = TRUE, outlier.size = 0.8) +
    geom_hline(yintercept = 0) +
    # separate by parameter
    facet_wrap(~ parameter+label, scales = "free", nrow = 1) + #, space = "free") +
    labs(x = NULL, y = "Relative error (%)", subtitle = "Distribution of relative error (RE; i.e., (true-estimated values)/true value*100)") +
    scale_x_discrete(labels = NULL, breaks = NULL)
}

p1pcod <- plot_sim(pcod_sim, "AI Pcod Simulation", fill_col = "goldenrod")
p1thrn <- plot_sim(thrn_sim, "GOA Thornyhead Simulation")

p2pcod <- plot_re(sim_re %>% filter(sp == "AI Pcod"), fill_col = "goldenrod")
p2thrn <- plot_re(sim_re %>% filter(sp == "GOA Thornyhead"))

cowplot::plot_grid(p1pcod +
                     labs(subtitle = "Distribution of parameters estimates\n(median=horizontal line, true value=point)"),
                   p2pcod +
                     labs(subtitle = "Distribution of relative error (RE;\ni.e., (true-estimated values)/true value*100)"),
                   ncol = 1)
ggsave(paste0("vignettes/ex4_sim_pcod.png"), width = 3.5, height = 6.5, units = "in", bg = "white")

cowplot::plot_grid(p1thrn, p2thrn, ncol = 1)
ggsave(paste0("vignettes/ex4_sim_thrn.png"), width = 11, height = 7, units = "in", bg = "white")

# mcmc ----

# function to (1) run models and (2) return posterior data frames and the traceplots
# input is model, number of iterations (samples), number of chains
mcmc_comp <- function(mod_name, it_numb, chain_numb) {

  # set up MCMC chain information
  it_num <- it_numb
  chain_num <- chain_numb
  # mod_name = pcod_mod; it_num = 4000; chain_num = 4

  # run model with laplace approximation
  mod_la <- tmbstan(obj = mod_name, chains = chain_num, init = mod_name$par, laplace = TRUE, iter = it_num)
  # run model without laplace approximation, i.e., all parameters fully estimated without assumptions of normality
  mod_mcmc <- tmbstan(obj = mod_name, chains = chain_num, init = mod_name$par, laplace = FALSE, iter = it_num)

  # posteriors as data frame
  post_la <- as.data.frame(mod_la); post_la$type <- ("la")
  post_mcmc <- as.data.frame(mod_mcmc); post_mcmc <- post_mcmc[, c(1:length(mod_name$par), dim(post_mcmc)[2])]; post_mcmc$type <- rep("mcmc")
  # informational things... this is for getting the posterior draws, i.e., to test traceplots and such
  post_la$chain <- rep(1:chain_num, each = it_num/2)
  post_la$iter_num <- rep(1:(it_num/2), chain_num)
  post_mcmc$chain <- rep(1:chain_num, each = it_num/2)
  post_mcmc$iter_num <- rep(1:(it_num/2), chain_num)
  post_draws <- rbind(post_la, post_mcmc)

  # get quantiles
  qv <- seq(from = 0, to = 1, by = 0.01)
  quant_dat <- data.frame(quant = NULL,
                          la = NULL,
                          mcmc = NULL,
                          par = NULL)

  for (i in 1:(dim(post_la)[2]-4)) { # post_la has type, chains, lp, iteration number column that don"t count

    tmp <- data.frame(quant = qv,
                      la = quantile(unlist(post_la[i]), probs = qv),
                      mcmc = quantile(unlist(post_mcmc[i]), probs = qv),
                      par = rep(paste0("V", i)))

    quant_dat <- rbind(quant_dat, tmp)

  }

  return(list(quant_dat, post_draws, mod_la, mod_mcmc))

}

# run models
pcod_comp <- mcmc_comp(mod_name = pcod_mod, it_numb = 1000, chain_numb = 3)
thrn_comp_log_tau <- mcmc_comp(mod_name = thrn_mod, it_numb = 5000, chain_numb = 3) # running slow bc of convergence issues

# parameter summary table - compare the marginal max likelihood estimates (MMLE)
# with the mcmc estimates. Rhat is the potential scale reduction factor on split
# chains (at convergence, Rhat=1)
m <- summary(pcod_comp[[3]])$summary[1, c(6,4,8,10), drop = FALSE]
psum <- data.frame(Stock = 'AI Pcod', Model = 'MCMC_Laplace', FixedEffect = row.names(m), m, row.names = NULL)
m <- summary(pcod_comp[[4]])$summary[1, c(6,4,8,10), drop = FALSE]
psum <- psum %>% bind_rows(data.frame(Stock = 'AI Pcod', Model = 'MCMC_withoutLaplace', FixedEffect = row.names(m), m, row.names = NULL))
m <- summary(thrn_comp_log_tau[[3]])$summary[1:6, c(6,4,8,10), drop = FALSE]
psum <- psum %>% bind_rows(data.frame(Stock = 'GOA Thornyhead', Model = 'MCMC_Laplace', FixedEffect = row.names(m), m, row.names = NULL))
m <- summary(thrn_comp_log_tau[[4]])$summary[1:6, c(6,4,8,10), drop = FALSE]
psum <- psum %>% bind_rows(data.frame(Stock = 'GOA Thornyhead', Model = 'MCMC_withoutLaplace', FixedEffect = row.names(m), m, row.names = NULL))
pars <- as.list(pcod_mod$sdrep, "Est")$log_PE
sd <- as.list(pcod_mod$sdrep, "Std")$log_PE
psum <- psum %>% bind_rows(data.frame(Stock = 'AI Pcod', Model = 'Base_MMLE_Laplace', FixedEffect = "log_PE",
                                      X50. = pars, X2.5. = pars-1.96*sd, X97.5. = pars+1.96*sd, Rhat = NA))
pars <- unlist(as.list(thrn_mod$sdrep, "Est"))
pars <- pars[names(pars)[grepl("log_PE|log_q|log_tau", names(pars))]]
sd <- unlist(as.list(thrn_mod$sdrep, "Std"))
sd <- sd[names(sd)[grepl("log_PE|log_q|log_tau", names(sd))]]
psum <- psum %>% bind_rows(data.frame(Stock = 'GOA Thornyhead', Model = 'Base_MMLE_Laplace', FixedEffect = rownames(summary(thrn_comp_log_tau[[3]])$summary[1:6,,drop=FALSE]),
                                      X50. = pars, sd = sd, row.names = NULL) %>%
                             mutate(X2.5. = X50.-1.96*sd,
                                    X97.5. = X50.+1.96*sd,
                                    Rhat = NA) %>%
                             select(-sd))
write.csv(psum, 'inst/example_data/ex4vignette_param_table.csv')
psum <- read.csv('inst/example_data/ex4vignette_param_table.csv')
psum <- psum %>%
  mutate(`Estimate (95% CI)` = ifelse(X50. > -1e3,
                                      # parameter table formatting: estimate (lci, uci)
                                      paste0(formatC(X50., format="f", big.mark=",", digits=2), " (",
                                             formatC(X2.5., format="f", big.mark=",", digits=2), ", ",
                                             formatC(X97.5., format="f", big.mark=",", digits=2), ")"),
                                      # if estimates are really small
                                      # (approaching 0 when exponentiated), use
                                      # scientific notation for the table
                                      paste0(formatC(X50., big.mark=",", digits=3), " (",
                                             formatC(X2.5., big.mark=",", digits=3), ", ",
                                             formatC(X97.5., big.mark=",", digits=3), ")")),
         Rhat = formatC(Rhat, format="f", digits=3)) %>%
           select(Stock, FixedEffect, Model, `Estimate (95% CI)`, Rhat) %>%
  arrange(Stock, FixedEffect, Model)

p1_mcmc_pcod <- rstan::traceplot(pcod_comp[[3]]) + ggtitle(label = "Traceplots to assess mixing across Markov chains", subtitle = "AI Pcod")
p1_mcmc_thrn <- rstan::traceplot(thrn_comp_log_tau[[3]], ncol = 3) + ggtitle(label = NULL, subtitle = "GOA Thornyhead")
cowplot::plot_grid(p1_mcmc_pcod, p1_mcmc_thrn, ncol = 1, rel_heights = c(1,1.5))
ggsave(paste0("vignettes/ex4_traceplots.png"), width = 8, height = 8, units = "in", bg = "white")

p2_mcmc_thrn <- bayesplot::mcmc_pairs(thrn_comp_log_tau[[3]],
                      off_diag_args = list(size = 0.8, alpha = 1/5),
                      pars = c("log_PE[1]", "log_PE[2]", "log_PE[3]",
                               "log_q", "log_tau_biomass", "log_tau_cpue"),
                      grid_args = list(top = "GOA Thornyhead: Pairwise correlation matrix of the posterior draws"))# +
p2_mcmc_thrn

# Pairwise correlation matrix of the posterior draws, with histograms of the
# univariate marginal distributions on the diagonal and a scatterplot of the
# bivariate distributions off the diagonal.

ggsave(plot = p2_mcmc_thrn, filename = paste0("vignettes/ex4_pairs_thrn.png"), width = 11, height = 8, units = "in", bg = "white")

# clean up data frame names by renaming things
pcod_qq <- pcod_comp[[1]]
pcod_qq$par_name <- rep("log_PE1"); pcod_qq$sp <- rep("AI Pcod")
thrn_qq <- thrn_comp_log_tau[[1]]
thrn_qq$par_name <- recode(thrn_qq$par,
                           V1 = "log_PE1",
                           V2 = "log_PE2",
                           V3 = "log_PE3",
                           V4 = "log_q",
                           V5 = "log_tau_biomass",
                           V6 = "log_tau_cpue")
thrn_qq$sp <- rep("GOA Thornyhead")
qq_dat <- rbind(pcod_qq, thrn_qq)

# plot
plot_laplace_mcmc <- function(qq_dat, plot_title) {
  ggplot(qq_dat, aes(mcmc, la)) +
    geom_abline(intercept = 0, slope = 1, col = "lightgray") +
    geom_point() +
    facet_wrap(~par_name, scales = "free", nrow = 2,  dir = "v") +
    labs(x = "MCMC", y = "MCMC with Laplace", title = plot_title) +
    theme(plot.title = element_text(size = 12),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 11),
          strip.text = element_text(size = 11))
}
p1 <- plot_laplace_mcmc(pcod_qq, "AI Pcod")
p1
ggsave(paste0("vignettes/ex4_mcmc_pcod_qq.png"), width = 3, height = 3, units = "in", bg = "white")
p2 <- plot_laplace_mcmc(thrn_qq, "GOA Thornyhead")
p2
ggsave(paste0("vignettes/ex4_mcmc_thrn_qq.png"), width = 7.5, height = 5, units = "in", bg = "white")

# cowplot::plot_grid(p1, p2, rel_widths = c(0.4,0.6), ncol = 2)
# ggsave(paste0("vignettes/ex4_mcmc_qq.png"), width = 11, height = 7, units = "in", bg = "white")

# OSAs -----

# Pcod
pcod_resid <- rema::get_osa_residuals(pcod_mod)
cowplot::plot_grid(pcod_resid$plots$qq +
                     ggtitle('OSA Residuals for AI Pcod') +
                     theme(legend.position = 'none'),
                   pcod_resid$plots$biomass_resids, ncol = 1)
ggsave(paste0("vignettes/ex4_aipcod_osa.png"), width = 5.5, height = 7, units = "in", bg = "white")

# Thorny
thrn_resid <- get_osa_residuals(thrn_mod)
p1 <-thrn_resid$plots$qq +
  ggtitle('OSA Residual QQ Plots for GOA Thornyhead') +
  theme(legend.position = 'none')
p2 <- thrn_resid$plots$biomass_qq + # biomass data
  ggtitle('By Biomass Survey Strata') +
  theme(legend.position = 'none')
p3 <- thrn_resid$plots$cpue_qq +
  ggtitle('By Longline Survey Strata') + # rpw data
  theme(legend.position = 'none')
cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1,2,1))
ggsave(paste0("vignettes/ex4_thrn_osaqq.png"), width = 7.5, height = 11, units = "in", bg = "white")

p1 <- thrn_resid$plots$biomass_resids + ggtitle('OSA Residuals for GOA Thornyhead by Biomass Survey Strata')
p2 <- thrn_resid$plots$cpue_resids + ggtitle('By Longline Survey Strata')
cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(2.5,1.5))
ggsave(paste0("vignettes/ex4_thrn_osa.png"), width = 7.5, height = 8, units = "in", bg = "white")

# likelihood profiling ----
# -   Is the model identifiable? Using a likelihood profiling approach, we will test if model parameters are identifiable, which means they have a solution (i.e., a minimum log-likelihood) and aren't highly correlated with another parameter
