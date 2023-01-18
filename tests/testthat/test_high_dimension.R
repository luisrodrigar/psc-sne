
n <- 100
x <- sphunif::r_unif_sph(n, 3, 3)
rho_list <- rep(0.5, n)
rho_list_diff <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), each = 10)

test_that("Checking that the calculus of high dimensional probabilities returns the expected value for both matrix and escalar versions.", {
  res_mat <- pscsne::high_dimension_mat(x, rho_list)
  res_esc <- pscsne::high_dimension(x, rho_list)
  expect_equal(res_esc, res_mat)
})

test_that("Checking that the calculus of high dimensional probabilities returns the same for both methods but in this case the rho_list contains different values.", {
  res_mat <- pscsne::high_dimension_mat(x, rho_list_diff)
  res_esc <- pscsne::high_dimension(x, rho_list_diff)
  expect_equal(res_esc, res_mat)
})
