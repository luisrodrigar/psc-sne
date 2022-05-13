library(rotasym)

simple_dspcauchy_sim <- function(s_ij, rho, d) {
  ((1 + rho^2 - 2 * rho * s_ij)^(-d))
}

q_i_not_j <- function(Y, i, j, rho, d) {
  n <- nrow(Y)
  qijs <- sapply(1:n, function(k){
    sapply(1:n, function(l) {
      return(simple_dspcauchy_ld(Y, l, k, rho, d))
    })
  }, simplify = 'array')
  diag(qijs) <- 0
  qijs[i, j] <- 0
  sum(qijs)
}

q_ij <- function(Y, s_ij, i, j, rho, d, has_deriv_num) {
  n <- nrow(Y)
  qij <- simple_dspcauchy_ld(Y, i, j, rho, d)
  if(has_deriv_num)
    qij <- simple_dspcauchy_sim(s_ij, rho, d)
  return(qij / Z(Y, s_ij, i, j, rho, d))
}

Z <- function(Y, s_ij, i, j, rho, d) {
  qij <- simple_dspcauchy_sim(s_ij, rho, d)
  q_inotj <- q_i_not_j(Y, i, j, rho, d)
  return(qij+q_inotj)
}

log_q_ij <- function(Y, s_ij, i, j, rho, d, has_deriv_num){
  log(q_ij(Y, s_ij, i, j, rho, d, has_deriv_num))
}

log_q_ij_prod_Z <- function(Y, s_ij, i, j, rho, d, has_deriv_num){
  log(q_ij(Y, s_ij, i, j, rho, d, has_deriv_num)*Z(Y, s_ij, i, j, rho, d))
}

log_Z <- function(Y, s_ij, i, j, rho, d){
  log(Z(Y, s_ij, i, j, rho, d))
}

simple_dspcauchy_yi <- function(Y, yi, j, rho, d) {
  ((1 + rho^2 - 2 * rho * t(yi) %*% Y[j,])^(-d))
}

simple_dspcauchy_yj <- function(Y, yi, i, rho, d) {
  ((1 + rho^2 - 2 * rho * t(Y[i, ]) %*% yi)^(-d))
}

simple_dspcauchy_ld <- function(Y, i, j, rho, d) {
  ((1 + rho^2 - 2 * rho * t(Y[i,]) %*% Y[j,])^(-d))
}

Z_yi <- function(Y, rho, d, yi, diff_i) {
  n <- nrow(Y)
  qijs <- sapply(1:n, function(k){
    sapply(1:n, function(l) {
      if(k==l) { 
        return(0)
      } else if(k == diff_i) {
        return(simple_dspcauchy_yi(Y, yi, l, rho, d))
      } else if(l == diff_i) {
        return(simple_dspcauchy_yj(Y, yi, k, rho, d))
      } else {
        return(simple_dspcauchy_ld(Y, k, l, rho, d))
      }
    })
  }, simplify = 'array')
  return(sum(qijs))
}

q_ij_yi <- function(Y, yi, i, j, rho, d, diff_i) {
  qij <- NULL
  if(diff_i == i) {
    qij <- simple_dspcauchy_yi(Y, yi, j, rho, d)
  } else if(diff_i == j) {
    qij <- simple_dspcauchy_yj(Y, yi, i, rho, d)
  } else {
    qij <- simple_dspcauchy_ld(Y, i, j, rho, d)
  }
  return(qij / Z_yi(Y, rho, d, yi, diff_i))
}

kl_cost_i <- function(Y, P, rho, d, yi, ii) {
  n <- nrow(Y)
  total <- sum(sapply(1:n, function(i){
    sapply((1:n), function(j) {
      if(i==j) {
        return(0)
      } else {
        (P[i,j]*log(P[i,j])-P[i,j] * log(q_ij_yi(Y, yi, i, j, rho, d, ii)))
      }
    })
  }))
  return(total)
}

gen_polysphere <- function(n, d, r) {
  p <- (d+1)
  polysphere <- array(NA, dim=c(n, p, r))
  for(k in seq_len(r)) {
    polysphere[,,k] <- r_unif_sphere(n, p)
  }
  polysphere
}

n <- 25
d <- 2
r <- 5
perplexity <- 15

X <- gen_polysphere(n, d, r)
rho_hat <- rho_optimize(X, perplexity)
P <- high_dimension_p(X, rho_hat, d)
P <- symmetric_probs(P)

rho <- 0.5
Y <- r_unif_sphere(n, d+1)
Q <- low_dimension_Q(Y, d, rho)

ii <- 1
yi <- Y[ii,]

test_that("Gradient Descent checking", {
  expect_equal(jacobian(kl_cost_i, Y=Y, P=P, rho=rho, d=d, x=yi, i=ii), kl_cost_i(Y, P, rho, d, yi, ii), tolerance=1e-4)
})
