
library(rotasym)
library(numDeriv)
library(mvtnorm)

simple_dspcauchy_sim <- function(s_ij, rho, d) {
  drop((1 + rho^2 - 2 * rho * s_ij)^(-d))
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
  drop((1 + rho^2 - 2 * rho * t(yi) %*% Y[j, ])^(-d))
}

simple_dspcauchy_yj <- function(Y, yi, i, rho, d) {
  drop((1 + rho^2 - 2 * rho * t(Y[i, ]) %*% yi)^(-d))
}

simple_dspcauchy_ld <- function(Y, i, j, rho, d) {
  drop((1 + rho^2 - 2 * rho * t(Y[i, ]) %*% Y[j, ])^(-d))
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
        (P[i, j] * log(P[i, j]) - P[i, j] *
           log(q_ij_yi(Y, yi, i, j, rho, d, ii)))
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

X <- sphunif::r_unif_sph(n, p + 1, r)
optim_hat <- rho_optim_bst(X, perplexity, num_cores = 2)
rho_hat <- optim_hat$rho_values
P <- optim_hat$P
P <- symmetric_probs(P)
Y <- sphunif::r_unif_sph(n, (d + 1))[, , 1]

ii <- sample(n, size = 1)
yi <- Y[ii, ]

test_that("Checking value with gradient approximation", {
  computation_grad_bar <- grad(func = function(x) {
    kl_div_obj_func(Y = Y, P = P, rho = rho, d = d, yi = x / sqrt(sum(x^2)),
                    ii = ii)
  }, x = yi, method = "simple", method.args = list(eps = 1e-6))
  expect_equal(computation_grad_bar, kl_divergence_grad(Y, ii, rho, d, P),
               tolerance = 1e-4)
})

test_that("Checking that the radial projection is done", {
  Y_non_in_the_sphere <- mvtnorm::rmvnorm(n, c(1, 2, 3), diag(3))
  Y_sphere <- radial_projection(Y_non_in_the_sphere)
  expect_equal(
    kl_divergence_grad(Y_non_in_the_sphere, ii, rho, d, P),
    kl_divergence_grad(Y_sphere, ii, rho, d, P)
  )
})

test_that("Checking that i is valid", {
  expect_error(kl_divergence_grad(Y, 120, rho, d, P))
})

test_that("Checking that the number of observations coincides with the number of rows in P", {
  Y_wrong_rows <- r_unif_sphere(2, p)
  expect_error(kl_divergence_grad(Y_wrong_rows, ii, rho, d, P))
})
