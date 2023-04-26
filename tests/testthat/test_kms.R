
n <- 50
mu <- rbind(c(1,0), c(0, 1))
kappa <- c(100, 75)
set.seed(123)
r_mvmf <- movMF::rmovMF(n, kappa * mu, c(0.6, 0.4))
h <- pscsne::bw_kms(r_mvmf)

set.seed(123)
r_mvmf_small <- movMF::rmovMF(n,
                              c(5, 10) * rbind(c(1, 0), c(1, 0)),
                              c(0.6, 0.4))
h_small <- pscsne::bw_kms(r_mvmf_small)

set.seed(123)
mu <- c(100, 150, 200) * rbind(c(1, 0), c(DirStats::to_cir(pi/2)),
                               (DirStats::to_cir(pi/4)))
r_mvmf_three <- movMF::rmovMF(n, mu, c(0.33, 0.33, 0.34))
h_three <- pscsne::bw_kms(r_mvmf_three)

test_that("Same value is returned with numDeriv and polykde::grad_hess_kde_polysph", {
  res_numDeriv <- kms_dir(r_mvmf, h = h, use_numderiv = FALSE)
  res_polykde <- kms_dir(r_mvmf, h = h, use_numderiv = TRUE)
  # Since numDeriv uses an approximation, the tolerance is defined to 1e-6
  expect_equal(res_numDeriv, res_polykde, tolerance = 1e-6)
})

test_that("From a mixture of two vMF with different mean and large concentration two groups are identified", {
  res <- kms_dir(r_mvmf, h = h)
  expect_equal(nrow(res$modes), 2)
})

test_that("From a mixture of two vMF with same mean and small concentration one group is identified", {
  res <- kms_dir(r_mvmf_small, h = h_small)
  expect_equal(nrow(res$modes), 1)
})

test_that("From a mixture of three vMF with different mean and high variable concentration parameter, then three groups are identified", {
  res <- kms_dir(r_mvmf_three, h = h_three)
  expect_equal(nrow(res$modes), 3)
})

test_that("From a mixture of three vMF with different mean and high variable concentration parameter, then two groups are identified when the number of initial clusters is passed as two", {
  res <- kms_dir(r_mvmf_three, h = h_three, num_init_clusters = 2)
  expect_equal(nrow(res$modes), 2)
})
