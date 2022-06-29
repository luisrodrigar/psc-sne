
p <- 2
r <- 3
n <- 3
i <- 1
j <- 2
rho <- 0.5
rho_list <- rep(rho, n)
k <- 1
set.seed(100469000)
x <- sphunif::r_unif_sph(n, p + 1, r)

test_that("HD: Checking method 'd_sph_cauchy' for a given i, j and k", {
  expected_value <- 2.58103
  expect_equal(d_sph_cauchy(x, i, j, rho, k, p),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("HD: Checking method 'd_psph_cauchy' for a given i and j", {
  expected_value <- 0.5543813
  expect_equal(d_psph_cauchy(x, i, j, rho),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("HD: Checking method 'd_i_psph_cauchy' for a given i", {
  expected_value <- 24.88156
  expect_equal(sum(d_i_psph_cauchy(x, i, rho_list)),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("HD: Checking method 'd_total_psph_cauchy' for each i-th observation", {
  expected_value <- c(24.8815586,  0.8615519, 24.6343479)
  expect_equal(d_total_psph_cauchy(x, rho_list),
    expected_value,
    tolerance = 1e-6, ignore_attr = TRUE
  )
})

test_that("HD: Checking method 'jcondi_psph', calculate the density function value for a given i to choose a given j observation as neighbor", {
  expected_value <- 0.02228081
  d_total_i_psph <- d_total_psph_cauchy(x, rho_list)[i]
  expect_equal(jcondi_psph(x, i, j, rho_list, d_total_i_psph),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("HD: Checking method 'prob_cond_i_psph', calculate probs for a given i to choose the j-th observation as neighbor", {
  expected_value <- c(0, 0.02228081, 0.97771919)
  expect_equal(prob_cond_i_psph(x, i, rho_list),
    expected_value,
    tolerance = 1e-6
  )
})

d <- 2
set.seed(100469000)
y <- sphunif::r_unif_sph(n, d + 1)[, , 1]

test_that("LD: Checking method 'd_sph_cauchy' for a given i and j observations", {
  expected_value <- 2.58103
  expect_equal(d_sph_cauchy(array(y, dim = c(nrow(y), ncol(y), 1)), i, j, rho, 1, d),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("LD: Checking method 'd_i_sph_cauchy' for a given i observation", {
  expected_value <- 18.06179
  expect_equal(d_i_sph_cauchy(y, i, rho),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("LD: Checking method 'jcondi_sph' for a given i and j observations", {
  expected_value <- 2.58103 / 18.06179
  expect_equal(jcondi_sph(y, i, j, rho),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("LD: Checking method 'prob_cond_i_sph' for a given i and j observations", {
  expected_value <- c(0, 2.58103 / 18.06179, 15.48076 / 18.06179)
  expect_equal(prob_cond_i_sph(y, i, rho),
    expected_value,
    tolerance = 1e-6, ignore_attr = TRUE
  )
})
