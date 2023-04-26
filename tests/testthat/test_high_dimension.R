
n <- 100
x <- sphunif::r_unif_sph(n, 3, 3)
rho_list <- rep(0.5, n)
rho_list_diff <- rep(1:10, each = 10)

test_that("High dimensional probabilities computed in matrix and scalar versions; common rho_list", {
  res_mat <- high_dim_mat(x, rho_list)
  res_esc <- high_dimension(x, rho_list)
  expect_equal(res_esc, res_mat)
})

test_that("High dimensional probabilities computed in matrix and scalar versions; rho_list with different values", {
  res_mat <- high_dim_mat(x, rho_list_diff)
  res_esc <- high_dimension(x, rho_list_diff)
  expect_equal(res_esc, res_mat)
})

test_that("high_dim_mat with cos_psh provided or not is the same", {
  n <- 7
  rho_list <- 1:n
  x <- sphunif::r_unif_sph(n = n, p = 3, M = 3)
  res1 <- high_dim_mat(x, rho_list)
  res2 <- high_dim_mat(x, rho_list, cos_psh = cosine_polysph(x))
  expect_equal(res1, res2)
})
