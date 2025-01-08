test_that("consistency in two-survey admb-re goasr_rwout.rep model results", {
  admb_re <- rema::read_admb_re(filename = 'testdata/goasr_rwout.rep')
  input <- rema::prepare_rema_input(admb_re = admb_re,
                                    multi_survey = 1)
  m <- rema::fit_rema(input)
  pred_biom <- rema::tidy_rema(m)$total_predicted_biomass
  # saveRDS(pred_biom, '_expect_goa_sr.rds')
  old_pred_biom <- readRDS('_expect_goa_sr.rds')
  expect_equal(pred_biom, old_pred_biom, tolerance = 0.1)
})
