
subtract_rows <- function(X, i) {
  n <- nrow(X)
  return(t(replicate(n - i, unlist(c(X[i, ])))) - X[-(1:i), ])
}

euclidean_distance <- function(X) {
  dist_val <- c()
  n <- nrow(X)
  for (i in seq_len(n - 1)) {
    sub_pwr_2 <- subtract_rows(X, i)^2
    sub_sqrt_sum <- sqrt(apply(sub_pwr_2, sum, MARGIN = 1))
    dist_val <- c(dist_val, rep(0, i), sub_sqrt_sum)
  }
  dist_val <- c(dist_val, rep(0, n))
  upper_triang <- matrix(dist_val, byrow = TRUE, nrow = n)
  return(t(upper_triang) + upper_triang)
}

index_except_i <- function(i, n) {
  index <- c(seq(1, i - 1), seq(i + 1, n))
  if (i == 1) {
    index <- 2:n
  } else if (i == n) {
    index <- 1:(n - 1)
  }
  return(index)
}

x_diff <- function(X) {
  sum_x <- apply(X^2, MARGIN = 1, FUN = sum)
  sum_x_m <- replicate(150, sum_x)
  cross_times_minus_2 <- -2 * (X %*% t(X))
  D <- t(cross_times_minus_2 + sum_x_m) + sum_x_m
  diag(D) <- 0
  return(D)
}

entropy_beta <- function(D_i, beta = 1) {
  P_i <- exp(-D_i * beta)
  sum_p_i <- sum(P_i)
  H_i <- log(sum_p_i) + (beta * sum(D_i * P_i) / sum_p_i)
  P_i <- P_i / sum_p_i
  return(list(entropy = H_i, probs = P_i))
}

binary_search <- function(h_diff, beta, i, beta_min, beta_max) {
  if (h_diff > 0) {
    beta_min <- beta[i]
    if (beta_max == -Inf || beta_max == Inf) {
      beta[i] <- beta[i] * 2
    } else {
      beta[i] <- (beta[i] + beta_max) / 2
    }
  } else {
    beta_max <- beta[i]
    if (beta_min == -Inf || beta_min == Inf) {
      beta[i] <- beta[i] / 2
    } else {
      beta[i] <- (beta[i] + beta_min) / 2
    }
  }
  return(list(beta = beta, min = beta_min, max = beta_max))
}

binary_search_optimization <- function(D_i, i, beta, h_star, prob_star, log_perp,
                                       tolerance = 1e-5) {
  beta_min <- -Inf
  beta_max <- Inf
  tries <- 0
  h_diff <- h_star - log_perp

  while (abs(h_diff) > tolerance && tries < 50) {
    beta_opt <- binary_search(h_diff, beta, i, beta_min, beta_max)
    beta <- beta_opt$beta
    beta_min <- beta_opt$min
    beta_max <- beta_opt$max

    res_loop <- entropy_beta(D_i, beta[i])
    h_star <- res_loop$entropy
    prob_star <- res_loop$probs

    h_diff <- h_star - log_perp
    tries <- tries + 1
  }
  return(list(probs = prob_star, beta = beta))
}
