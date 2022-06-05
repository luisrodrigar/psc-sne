library(rotasym)
library(numDeriv)
library(mvtnorm)

simple_dspcauchy_sim <- function(s_ij, rho, d) {
  ((1 + rho^2 - 2 * rho * s_ij)^(-d))
}

q_i_not_j <- function(Y, i, j, rho, d) {
  n <- nrow(Y)
  qijs <- sapply(1:n, function(k) {
    sapply(1:n, function(l) {
      return(simple_dspcauchy_ld(Y, l, k, rho, d))
    })
  }, simplify = "array")
  diag(qijs) <- 0
  qijs[i, j] <- 0
  sum(qijs)
}

q_ij <- function(Y, s_ij, i, j, rho, d, has_deriv_num) {
  n <- nrow(Y)
  qij <- simple_dspcauchy_ld(Y, i, j, rho, d)
  if (has_deriv_num) {
    qij <- simple_dspcauchy_sim(s_ij, rho, d)
  }
  return(qij / Z(Y, s_ij, i, j, rho, d))
}

Z <- function(Y, s_ij, i, j, rho, d) {
  qij <- simple_dspcauchy_sim(s_ij, rho, d)
  q_inotj <- q_i_not_j(Y, i, j, rho, d)
  return(qij + q_inotj)
}

log_q_ij <- function(Y, s_ij, i, j, rho, d, has_deriv_num) {
  log(q_ij(Y, s_ij, i, j, rho, d, has_deriv_num))
}

log_q_ij_prod_Z <- function(Y, s_ij, i, j, rho, d, has_deriv_num) {
  log(q_ij(Y, s_ij, i, j, rho, d, has_deriv_num) * Z(Y, s_ij, i, j, rho, d))
}

log_Z <- function(Y, s_ij, i, j, rho, d) {
  log(Z(Y, s_ij, i, j, rho, d))
}

simple_dspcauchy_yi <- function(Y, yi, j, rho, d) {
  ((1 + rho^2 - 2 * rho * t(yi) %*% Y[j, ])^(-d))
}

simple_dspcauchy_yj <- function(Y, yi, i, rho, d) {
  ((1 + rho^2 - 2 * rho * t(Y[i, ]) %*% yi)^(-d))
}

simple_dspcauchy_ld <- function(Y, i, j, rho, d) {
  ((1 + rho^2 - 2 * rho * t(Y[i, ]) %*% Y[j, ])^(-d))
}

Z_yi <- function(Y, rho, d, yi, diff_i) {
  n <- nrow(Y)
  qijs <- sapply(1:n, function(k) {
    sapply(1:n, function(l) {
      if (k == l) {
        return(0)
      } else if (k == diff_i) {
        return(simple_dspcauchy_yi(Y, yi, l, rho, d))
      } else if (l == diff_i) {
        return(simple_dspcauchy_yj(Y, yi, k, rho, d))
      } else {
        return(simple_dspcauchy_ld(Y, k, l, rho, d))
      }
    })
  }, simplify = "array")
  return(sum(qijs))
}

q_ij_yi <- function(Y, yi, i, j, rho, d, diff_i) {
  qij <- NULL
  if (diff_i == i) {
    qij <- simple_dspcauchy_yi(Y, yi, j, rho, d)
  } else if (diff_i == j) {
    qij <- simple_dspcauchy_yj(Y, yi, i, rho, d)
  } else {
    qij <- simple_dspcauchy_ld(Y, i, j, rho, d)
  }
  return(qij / Z_yi(Y, rho, d, yi, diff_i))
}

kl_div_obj_func <- function(Y, P, rho, d, yi, ii) {
  Y <- radial_projection(Y)
  n <- nrow(Y)
  total <- sum(sapply(1:n, function(i) {
    sapply((1:n), function(j) {
      if (i == j) {
        return(0)
      } else {
        (P[i, j] * log(P[i, j]) - P[i, j] * log(q_ij_yi(Y, yi, i, j, rho, d, ii)))
      }
    })
  }))
  return(total)
}

n <- 25
p <- 2
d <- 2
r <- 5
rho <- 0.5
perplexity <- 15

X <- gen_polysphere(n, p, r)
rho_hat <- rho_optim_bst(X, perplexity)$rho_values
P <- high_dimension(X, rho_hat)
P <- symmetric_probs(P)
Y <- r_unif_sphere(n, (d + 1))

ii <- 1
yi <- Y[ii, ]

test_that("Checking value with gradient approximation", {
  expect_equal(jacobian(kl_div_obj_func, Y = Y, P = P, rho = rho, d = d, x = yi, i = ii),
    kl_divergence_grad(Y, ii, rho, d, P),
    tolerance = 1e-6, ignore_attr = TRUE
  )
})

test_that("Checking that the radial projetion is done", {
  Y_non_in_the_sphere <- rmvnorm(n, c(1, 2, 3), diag(3))
  Y_sphere <- radial_projection(Y_non_in_the_sphere)
  expect_equal(
    kl_divergence_grad(Y_non_in_the_sphere, ii, rho, d, P),
    kl_divergence_grad(Y_sphere, ii, rho, d, P)
  )
})

test_that("Checking that i is valid", {
  expect_error(kl_divergence_grad(Y, 120, rho, d, P))
})

test_that("Checking that d is valid", {
  expect_error(kl_divergence_grad(Y, ii, rho, 0, P))
})

test_that("Checking that rho is valid", {
  expect_error(kl_divergence_grad(Y, ii, 1, d, P))
  expect_error(kl_divergence_grad(Y, ii, -1, d, P))
})

test_that("Checking that the number of observation is according with the number of rows in P", {
  Y_wrong_rows <- r_unif_sphere(2, p)
  expect_error(kl_divergence_grad(Y_wrong_rows, ii, rho, d, P))
})

test_that("Checking that the dimension of Y is according with the given d", {
  Y_wrong_rows <- r_unif_sphere(n, 1)
  expect_error(kl_divergence_grad(Y_wrong_rows, ii, rho, d, P))
  Y_wrong_rows <- r_unif_sphere(n, 2)
  expect_error(kl_divergence_grad(Y_wrong_rows, ii, rho, d, P))
})
