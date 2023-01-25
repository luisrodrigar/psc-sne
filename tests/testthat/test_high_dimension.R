
n <- 100
x <- sphunif::r_unif_sph(n, 3, 3)
rho_list <- rep(0.5, n)
rho_list_diff <- rep(seq(1, 10), each = 10)

test_that("High dimensional probabilities computed in matrix and scalar versions; common rho_list", {
  res_mat <- pscsne::high_dim_mat(x, rho_list)
  res_esc <- pscsne::high_dimension(x, rho_list)
  expect_equal(res_esc, res_mat)
})

test_that("High dimensional probabilities computed in matrix and scalar versions; rho_list with different values", {
  res_mat <- pscsne::high_dim_mat(x, rho_list_diff)
  res_esc <- pscsne::high_dimension(x, rho_list_diff)
  expect_equal(res_esc, res_mat)
})
