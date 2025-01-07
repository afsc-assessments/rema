test_that("consistency in extra_cv using two-survey admb-re goasr_rwout.rep model results", {
  admb_re <- rema::read_admb_re(filename = 'testdata/goasr_rwout.rep')
  input <- rema::prepare_rema_input(admb_re = admb_re,
                                    multi_survey = 1,
                                    extra_biomass_cv = list(assumption = 'extra_cv'))
  m <- rema::fit_rema(input)
  params <- rema::tidy_rema(m)$parameter_estimates
  # saveRDS(params, '_expect_extracv_goa_sr.rds')
  old_params <- readRDS('_expect_extracv_goa_sr.rds')
  testthat::expect_equal(params, old_params)
})
