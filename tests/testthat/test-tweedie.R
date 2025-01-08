test_that("consistency in tweedie results using ebs shelf", {
  biomass_dat <- read.csv('testdata/ebsshelf_orox.csv')
  input <- rema::prepare_rema_input(biomass_dat = biomass_dat,
                                    zeros = list(assumption = 'tweedie'))
  m <- rema::fit_rema(input)
  pred_biom <- rema::tidy_rema(m)$total_predicted_biomass
  # saveRDS(pred_biom, '_expect_tweedie_ebsshelf_orox.rds')
  old_pred_biom <- readRDS('_expect_tweedie_ebsshelf_orox.rds')
  expect_equal(pred_biom, old_pred_biom, tolerance = 0.1)
})
