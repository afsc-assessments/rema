test_that("consistency in univariate admb-re ai_shortraker_rwout.rep model results", {
  admb_re <- rema::read_admb_re(filename = 'testdata/ai_shortraker_rwout.rep')
  input <- rema::prepare_rema_input(admb_re = admb_re)
  m <- rema::fit_rema(input)
  pred_biom <- rema::tidy_rema(m)$total_predicted_biomass
  # saveRDS(pred_biom, '_expect_ai_shortraker.rds')
  old_pred_biom <- readRDS('_expect_ai_shortraker.rds')
  expect_equal(pred_biom, old_pred_biom)
  })
