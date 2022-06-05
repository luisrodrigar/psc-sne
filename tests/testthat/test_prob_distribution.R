
p <- 2
r <- 3
n <- 3
i <- 1
j <- 2
rho <- 0.5
rho_list <- rep(rho, n)
k <- 1
set.seed(100469000)
x <- gen_polysphere(n, p, r)

test_that("HD: Checking method 'simple_dspcauchy_hd' for a given i, j and k", {
  expected_value <- 0.241892
  expect_equal(simple_dspcauchy_hd(x, i, j, rho, k, p),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("HD: Checking method 'P_ij_psc' for a given i and j", {
  expected_value <- 0.03305185
  expect_equal(P_ij_psc(x, i, j, rho),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("HD: Checking method 'P_i_psc' for a given i", {
  expected_value <- 0.07566693
  expect_equal(sum(P_i_psc(x, i, rho_list)),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("HD: Checking method 'P_total_psc' for each i-th observation", {
  expected_value <- c(0.07566693, 0.9713986, 0.9809618)
  expect_equal(P_total_psc(x, rho_list),
    expected_value,
    tolerance = 1e-6, ignore_attr = TRUE
  )
})

test_that("HD: Checking method 'jcondi_psc', calculate prob for a given i to choose a given j observation as neighbor", {
  expected_value <- 0.03305185 / 0.07566693
  expect_equal(jcondi_psc(x, i, j, rho_list, 0.07566693),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("HD: Checking method 'psc_cond_given_i', calculate probs for a given i to choose the j-th observation as neighbor", {
  expected_value <- c(0, 0.03305185 / 0.07566693, 0.04261508 / 0.07566693)
  expect_equal(psc_cond_given_i(x, i, rho_list, 0.07566693),
    expected_value,
    tolerance = 1e-6
  )
})

d <- 2
set.seed(100469000)
y <- sphunif::r_unif_sph(n, d + 1)[, , 1]

test_that("LD: Checking method 'simple_dspcauchy_ld' for a given i and j observations", {
  expected_value <- 2.58103
  expect_equal(simple_dspcauchy_ld(y, i, j, rho, d),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("LD: Checking method 'Q_i_sc' for a given i observation", {
  expected_value <- 18.06179
  expect_equal(Q_i_sc(y, i, rho),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("LD: Checking method 'jcondi_sc' for a given i and j observations", {
  expected_value <- 2.58103 / 18.06179
  expect_equal(jcondi_sc(y, i, j, rho),
    expected_value,
    tolerance = 1e-6
  )
})

test_that("LD: Checking method 'jcondi_sc' for a given i and j observations", {
  expected_value <- c(0, 2.58103 / 18.06179, 15.48076 / 18.06179)
  expect_equal(sc_cond_given_i(y, i, rho),
    expected_value,
    tolerance = 1e-6, ignore_attr = TRUE
  )
})
