test_that("consistency in multivariate admb-re bsai_sst_rwout.rep model results", {
  admb_re <- rema::read_admb_re(filename = 'testdata/bsaisst_rwout.rep')
  input <- rema::prepare_rema_input(admb_re = admb_re)
  m <- rema::fit_rema(input)
  pred_biom <- rema::tidy_rema(m)$total_predicted_biomass
  # saveRDS(pred_biom, '_expect_bsai_sst.rds')
  old_pred_biom <- readRDS('_expect_bsai_sst.rds')
  expect_equal(pred_biom, old_pred_biom, tolerance = 0.1)
})
