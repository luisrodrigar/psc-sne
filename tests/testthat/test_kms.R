
mu <- rbind(c(1,0), c(0,1))
kappa <- c(100, 75)
set.seed(123)
r_mvmf <- movMF::rmovMF(10, kappa * mu, c(0.6, 0.4))
h <- pscsne::bw_kms(r_mvmf)

test_that("Same value is returned with numDeriv and polykde::grad_hess_kde_polysph", {
  res_numDeriv <- kms_dir(r_mvmf, h = h, is_numDeriv = FALSE)
  res_polykde <- kms_dir(r_mvmf, h = h, is_numDeriv = TRUE)
  # Since numDeriv uses an approximation, the tolerance is defined to 1e-6
  expect_equal(res_numDeriv, res_polykde, tolerance = 1e-6)
})

test_that("From a mixture of two vMF two groups are identified", {
  res <- kms_dir(r_mvmf, h = h)
  expect_equal(nrow(res$modes), 2)
})
