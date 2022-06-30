###########################################
##      spherical Cauchy Perplexity      ##
###########################################

#' @title Perplexity of the \eqn{i}-th observation (scalar version)
#'
#' @description Calculate the perplexity of the \eqn{i}-th observation for a given a rho parameter.
#'
#' @inheritParams high_dimension
#' @param i corresponds to the \eqn{i}-th observation for the perplexity is calculated.
#' @param rho concentration parameter must be in [0, 1).
#' @return Perplexity and probabilities of the \eqn{i}-th observation for all the remainder observations.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(25, 3, 3)
#' to_perplexity_P(x, 1, 0.5)
#' to_perplexity_P(x, 4, 1 - 1e-4)
to_perplexity_P <- function(x, i, rho) {
  if (i < 1 || i > nrow(x)) {
    stop("i not valid, must be in [1, nrow(x)]")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("rho must be an scalar")
  }
  # Calculate the total marginal probability for the i-th observation
  d_total_i_psph <- sum(d_i_psph_cauchy(x = x, i = i, rho_list = rep(rho, nrow(x))))
  # Calculate the probability of the j-th given the i-th observation
  Pjcondi <- prob_cond_i_psph(
    x = x, i = i, rho_list = rep(rho, nrow(x)),
    d_total_i_psph_cauchy = d_total_i_psph
  )
  # Entropy formula for the j-th observation given the i-th one
  entropy <- function(j) {
    return(Pjcondi[j] * log2(Pjcondi[j]))
  }
  # Calculate the perplexity
  perp <- 2^(-sum(sapply(seq_len(nrow(x))[-i], entropy)))
  return(list(perp = perp, Pjcondi = Pjcondi))
}

#' @title Perplexity of the \eqn{i}-th observation (semi-scalar version)
#'
#' @description Calculate the perplexity of the \eqn{i}-th observation for a given a rho \eqn{\rho} concentration parameter.
#'
#' @inheritParams high_dimension
#' @inheritParams to_perplexity_P
#' @param cos_sim_ps cosine similarities matrix of each sphere \eqn{\mathcal{S}^p} (optional parameter, default value \code{NULL}).
#' @return Perplexity of the \eqn{i}-th observation for all the remainder observations.
to_perplexity <- function(x, i, rho, cos_sim_ps = NULL) {
  if (i < 1 || i > nrow(x)) {
    stop("i not valid, must be in [1, nrow(x)]")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("rho must be an scalar")
  }
  if (!is.null(cos_sim_ps) && length(dim(cosine_polysph(x))) != 3) {
    stop("cos_sim_ps has to be an array of size c(n, p + 1, r), from (S^p)^r")
  }
  # Calculate the cosine similarities of 'x' if 'cos_sim_ps' param is null
  if (is.null(cos_sim_ps)) {
    cos_sim_ps <- cosine_polysph(x)
  }
  # Conditional probability of the high-dimension of x
  Pjcondi <- high_dimension(x, rep(rho, nrow(x)), cos_sim_ps)
  # Entropy formula for the j-th observation given the i-th one
  entropy <- function(j) {
    (Pjcondi[i, j] * log2(Pjcondi[i, j]))
  }
  return(2^(-sum(sapply(seq_len(nrow(x))[-i], entropy))))
}

# check <- function(l) max(sapply(l, function(y) max(abs(l[[1]] - y)))) < 1e-7
# microbenchmark::microbenchmark(
#   to_perplexity_P(x, 1, 0.5),
#   to_perplexity(x, 1, 0.5),
#   check = check
# )
#
# Unit: milliseconds
# expr       min        lq      mean    median        uq       max neval
# to_perplexity_P(x, 1, 0.5)  16.69372  17.88449  24.63309  18.65393  27.29528  234.1621   100
# to_perplexity(x, 1, 0.5) 517.03870 583.15446 625.63860 613.75799 638.68714 1649.2364   100

#' @title Perplexity matrix (scalar version)
#'
#' @description Calculate the perplexity of each observations for a given rho parameters list.
#'
#' @inheritParams high_dimension
#' @return Perplexity of the \eqn{i}-th observation for all the remainder observations.
to_perp_scalar <- function(x, rho_list) {
  if (length(dim(x)) != 3) {
    stop("x must be an array of dimension c(n, p + 1, r), from (S^p)^r")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("rho_list must be a vector")
  }
  # Call to_perplexity_P (scalar way) for all the observations
  return(sapply(
    X = 1:nrow(x),
    FUN = function(i, x, rho_list) to_perplexity_P(x, i, rho_list[i])$perp,
    x = x, rho_list = rho_list
  ))
}

#' @title Perplexity matrix (matrix version)
#'
#' @description Calculate the perplexity of each observations for a given \eqn{\mathbf{\rho}} parameters list. Matrix version algorithm.
#'
#' @inheritParams high_dimension
#' @inheritParams to_perplexity
#' @return Perplexity of the \eqn{i}-th observation for all the remainder observations.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(25, 3, 3)
#' to_perp(x, rep(0.5, 25))
#' to_perp(x, rep(1 - 1e-4, 25), cosine_polysph(x))
to_perp <- function(x, rho_list, cos_sim_ps = NULL) {
  if (length(dim(x)) != 3) {
    stop("x must be an array of dimension c(n, p + 1, r), from (S^p)^r")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("rho_list must be a vector")
  }
  if (!is.null(cos_sim_ps) && length(dim(cos_sim_ps)) != 3) {
    stop("cos_sim_ps must be an array of dimension c(n, n, r) where (S^p)^r")
  }
  # Calculate the cosine similarities of 'x' if 'cos_sim_ps' param is null
  if (is.null(cos_sim_ps)) {
    cos_sim_ps <- cosine_polysph(x)
  }
  # Calculate the high-dimension probabilities
  P <- high_dimension(x, rho_list, cos_sim_ps)
  # Apply the formula of the entropy (matrix way)
  op <- P * log2(P)
  # Set the diagonal elements to zero
  diag(op) <- 0
  # Return the perplexity list result after applying the perplexity formula (matrix way)
  return(2^(-rowSums(op)))
}

# microbenchmark::microbenchmark(
#   to_perplexity_P(x, 1, 0.5),
#   to_perp(x, rep(.5, nrow(x)), cosine_polysph(x))
# )
#
# Unit: milliseconds
# expr      min       lq      mean    median        uq      max neval
# to_perplexity_P(x, 1, 0.5) 16.52357  18.1450  26.02556  19.03716  32.80008 175.8253   100
# to_perp(x, rep(0.5, nrow(x)), cos_sim_sph) 96.76836 123.6451 152.28035 142.41592 149.62526 381.0817   100

#####################################################
##  Optimize rho list based on a fixed perplexity  ##
#####################################################

### scalar calculus

# Time difference of 10.21872 mins

#' @title Serial optimization \eqn{\rho} concentration parameters (scalar version)
#'
#' @description Calculate the rho list values based on a fixed perplexity and a given data in \eqn{(\mathcal{S}^p)^r}.
#' Optimize the value using a method L-BFGS-B, setting the boundaries from 0 to 0.9999.
#' Each value is calculated serially. It prints the time consumption.
#'
#' @inheritParams high_dimension
#' @param perplexity a fixed value (between 5 and 100) to optimize the rho parameters.
#' @return Rho list (\eqn{\boldsymbol{\rho}}) with the values optimized for the given perplexity.
rho_optim_serial <- function(x, perplexity) {
  # Sample size
  n <- nrow(x)
  # Time start
  start_time <- Sys.time()
  # For each observation, calculate the optimal value based on the fixed perplexity
  rho_opt <- sapply(1:n, function(i) {
    stats::optim(
      par = 0.5,
      fn = function(rho) {
        # Objective function: (perplexity-fixed_perplexity)^2
        res <- (to_perplexity_P(x, i, rho)$perp - perplexity)^2
        ifelse(is.finite(res), res, 1e6)
      },
      method = "L-BFGS-B", lower = 0, upper = 1 - 1e-4
    )$par
  })
  # Time end
  end_time <- Sys.time()
  # Print time difference
  print(end_time - start_time)
  # Simplifying the result
  rho_opt <- simplify2array(rho_opt)
  return(rho_opt)
}

# Time difference of 18.12201 mins

#' @title Concurrent optimization of the \eqn{\rho} concentration parameters (scalar version)
#'
#' @description Calculate the rho list values based on a fixed perplexity and a given data in \eqn{(\mathcal{S}^p)^r}.
#' Optimize the value using the method L-BFGS-B, setting the boundaries within [0, 1).
#' The limit is considered as \code{1 - 1e-4}. Each value is calculated concurrently. At the end, the time consumption is shown.
#'
#' @inheritParams high_dimension
#' @inheritParams rho_optim_serial
#' @param num_cores number of cores to execute the code concurrently. This value must be below the total number of the CPU has available.
#' @return Rho list (\eqn{\boldsymbol{\rho}}) with the values optimized for the given perplexity.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(20, 3, 4)
#' rho_optim_par(x, 22, 2)
#' rho_optim_par(x, 30, 2)
rho_optim_par <- function(x, perplexity, num_cores = parallel::detectCores() - 1) {
  # Sample size
  n <- nrow(x)
  # Setting up the characteristics of the parallelization
  cl <- clusterFactory(num_cores)
  # Time start
  start_time <- Sys.time()
  # For each observation, calculate concurrently the optimal value based on the fixed perplexity
  rho_opt <- parallel::parLapply(cl, 1:n, function(i) {
    stats::optim(
      par = 0.5,
      fn = function(rho) {
        # Objective function: (perplexity-fixed_perplexity)^2
        res <- (to_perplexity_P(x, i, rho)$perp - perplexity)^2
        ifelse(is.finite(res), res, 1e6)
      },
      method = "L-BFGS-B", lower = 0, upper = 1 - 1e-4
    )$par
  })
  # Time end
  end_time <- Sys.time()
  # Print time difference
  print(end_time - start_time)
  # Simplifying the result
  rho_opt <- simplify2array(rho_opt)
  # Stop the clusters
  parallel::stopCluster(cl)
  return(rho_opt)
}

# > system.time(rho_optim_par(spokes, 22))
# Time difference of 29.04886 mins
# user   system  elapsed
# 3.812    5.742 1743.000

#' @title Concurrent optim to calculate \eqn{\rho} concentration parameters (scalar version)
#'
#' @description Calculate the rho list values based on a fixed perplexity and a given data in \eqn{(\mathcal{S}^p)^r}.
#' Optimize the value using the concurrently method L-BFGS-B (optimParallel).
#' The boundaries are set from 0 to 0.9999. It prints the time consumption.
#'
#' @inheritParams high_dimension
#' @inheritParams rho_optim_serial
#' @inheritParams rho_optim_par
#' @return Rho list (\eqn{\boldsymbol{\rho}}) with the values optimized for the given perplexity.
rho_optimParallel <- function(x, perplexity, num_cores = parallel::detectCores() - 1) {
  # Sample size
  n <- nrow(x)
  # Setting up the characteristics of the parallelization
  cl <- clusterFactory(num_cores)
  # Time start
  start_time <- Sys.time()
  # For each observation, calculate the optimal value using a concurrently optimization method
  rho_opt <- sapply(1:n, FUN = function(i) {
    optimParallel::optimParallel(
      par = 0.5,
      fn = function(rho) {
        # Objective function: (perplexity-fixed_perplexity)^2
        res <- to_perplexity_P(x, i, rho)$perp - perplexity
        ifelse(is.finite(res), res, 1e6)
      },
      lower = 0, upper = 1 - 1e-4,
      parallel = list(cl = cl, forward = FALSE)
    )$par
  })
  # Time end
  end_time <- Sys.time()
  # Print time difference
  print(end_time - start_time)
  # Stop the clusters
  parallel::stopCluster(cl)
  return(rho_opt)
}

# system.time(rho_optimParallel(spokes, 22))
# Time difference of 1.281991 hours
# user   system  elapsed
# 23.811   15.432 4615.212

## matricidal calculus, cost function
## Time difference of 14.22428 mins
## 50 rows 25 spheres

#' @title Serial optimization of the \eqn{\rho} concentration parameters (matrix version)
#'
#' @description Calculate the rho list values based on a fixed perplexity and a given data in \eqn{(\mathcal{S}^p)^r}.
#' Optimize the value using the method L-BFGS-B. The boundaries are set from 0 to 0.9999. It prints the time consumption.
#'
#' @inheritParams high_dimension
#' @inheritParams rho_optim_serial
#' @return Rho list (\eqn{\boldsymbol{\rho}}) with the values optimized for the given perplexity.
rho_optim_ineff <- function(x, perplexity) {
  # Sample size
  n <- nrow(x)
  # Time start
  start_time <- Sys.time()
  # Calculate the cosine similarities of (S^p)^r
  cosine_polysph <- cosine_polysph(x)
  # For each observation, calculate the optimal value matrix way
  rho_opt <- sapply(1:n, function(i) {
    stats::optim(
      par = rep(0.5, n),
      fn = function(rho) {
        # Objective function |perplexity-fixed_perplexity|_2
        dif <- to_perp(x, rho, cosine_polysph) - perplexity
        res <- t(dif) %*% dif
        ifelse(is.finite(res), res, 1e6)
      },
      method = "L-BFGS-B", lower = rep(0, n), upper = rep(1 - 1e-4, n)
    )$par
  })
  # Time end
  end_time <- Sys.time()
  # Print time difference
  print(end_time - start_time)
  return(rho_opt)
}

## efficient ways to calculate the perplexity

## parallel and matrix calculus
## Time difference of 14.56007 mins
## 50 rows 25 spheres

#' @title Concurrent optimization of the \eqn{\rho} concentration parameters (matrix version)
#'
#' @description Calculate the rho list values based on a fixed perplexity and a given data in \eqn{(\mathcal{S}^p)^r}.
#' Optimize the value using the method L-BFGS-B.
#' The boundaries are set from 0 to 0.9999. It prints the time consumption.
#'
#' @inheritParams high_dimension
#' @inheritParams rho_optim_serial
#' @inheritParams rho_optim_par
#' @return Rho list (\eqn{\boldsymbol{\rho}}) with the values optimized for the given perplexity.
rho_optimize_1 <- function(x, perplexity, num_cores = parallel::detectCores() - 1) {
  # Sample size
  n <- nrow(x)
  # Setting up the characteristics of the parallelization
  cl <- parallel::makeForkCluster(num_cores)
  # Time start
  start_time <- Sys.time()
  # Calculate the cosine similarities of (S^p)^r
  cosine_polysph <- cosine_polysph(x)
  # For each observation, calculate concurrently the optimal value based on the fixed perplexity
  rho_opt <- parallel::parLapply(cl, 1:n, function(i) {
    stats::optim(
      par = 0.5,
      fn = function(rho) {
        # Objective function: (perplexity-fixed_perplexity)^2
        res <- (to_perplexity(x = x, i = i, rho = rho, cosine_polysph) - perplexity)^2
        ifelse(is.finite(res), res, 1e6)
      },
      method = "L-BFGS-B", lower = 0, upper = 1 - 1e-4
    )$par
  })
  # Time end
  end_time <- Sys.time()
  # Print time difference
  print(end_time - start_time)
  # Simplifying the result
  rho_opt <- simplify2array(rho_opt)
  # Stop the clusters
  parallel::stopCluster(cl) # Time difference of 5.836471 mins
  return(rho_opt)
}

# check <- function(l) max(sapply(l, function(y) max(abs(l[[1]] - y)))) < 1e-7
# microbenchmark::microbenchmark(
#  rho_optimize(xxx, 15),
#  rho_optimize_2(xxx, 15),
#  check = check
# )

# Unit: milliseconds
# expr        min         lq       mean     median         uq        max neval
# rho_optimize(xxx, 15)   64.95212   74.77625   89.47905   84.41412   95.07279   200.1545   100
# rho_optimize_2(xxx, 15) 7582.09367 7729.17412 8035.37894 7880.25502 8186.76657 11580.5412   100

# check <- function(l) max(sapply(l, function(y) max(abs(l[[1]] - y)))) < 1e-7
# microbenchmark::microbenchmark(
#  rho_optimize(xxx, 15),
#  rho_optim_ineff(xxx, 15),
#  check = check
# )

#' @title  Binary search tree algorithm
#'
#' @description Calculate the rho value based on a perplexity difference
#' and current values of rho, rho min and rho max.
#'
#' @param perp_diff difference between the value obtained and the fixed perplexity
#' @param rho the current value optimized
#' @param rho_min min value for the current rho
#' @param rho_max max value for the current rho
#' @return Rho concentration parameter found in this step, currently min and max rho
bin_search <- function(perp_diff, rho, rho_min, rho_max) {
  if (perp_diff > 0 || is.na(perp_diff)) {
    rho_min <- rho
  } else {
    rho_max <- rho
  }
  rho <- (rho_min + rho_max) / 2
  return(list(rho = rho, min = rho_min, max = rho_max))
}

#' @title Binary search rho optimization for the \eqn{i}-th observation
#'
#' @description  Calculate the rho based on a fixed perplexity and a given data in \eqn{(\mathcal{S}^p)^r}.
#' The boundaries are set from 0 to 0.9999.
#'
#' @inheritParams high_dimension
#' @param i the \eqn{i}-th observation.
#' @param perp_fixed a fixed value used to optimized the rho values.
#' @param tolerance whether the difference between previous and current results is below this value (optional, default value \code{1e-3}).
#' @param rho parameter which determines the concentration of the spherical Cauchy distribution (optional, default value \code{0.5}).
#' @param max_tries number of maximum tries for each value of the rho list.
#' @return Rho concentration parameter and conditional probabilities calculated for the \eqn{i}-th observation.
rho_optim_i_bst <- function(x, i, perp_fixed, tolerance = 1e-3, rho = 0.5,
                            max_tries = 20) {
  # Min value of rho
  rho_min <- 0
  # Max value of rho
  rho_max <- 1 - 1e-4
  # Number of tries
  tries <- 0
  # Difference of perplexity, initially set to NA
  perp_diff <- NA

  # Calculate perplexity based on the default value of rho
  res <- to_perplexity_P(x, i, rho)
  # Get the perplexity value
  perp_star <- res$perp
  # Get the conditional probabilities given i
  prob_star <- res$Pjcondi

  # Set the difference to - infinity if perplexity is NA, calculate difference otherwise
  if (is.na(perp_star)) {
    perp_diff <- -Inf
  } else {
    perp_diff <- perp_star - perp_fixed
  }

  # While difference of perplexity is NA or greater than tolerance
  # and the number of tries less than the maximum, do the following
  while ((is.na(perp_diff) || abs(perp_diff) > tolerance) && tries <= max_tries) {
    # Find by means of the binary search optimal rho
    rho_opt <- bin_search(perp_diff, rho, rho_min, rho_max)
    # Get the current, min and max rho values
    rho <- rho_opt$rho
    rho_min <- rho_opt$min
    rho_max <- rho_opt$max

    # Calculate the perplexity and the P_j|i for the values of i-th iteration
    res_loop <- to_perplexity_P(x, i, rho)
    perp_star <- res_loop$perp
    prob_star <- res_loop$Pjcondi

    # Again, calculate the difference of perplexity as it was done outside the loop
    if (is.na(perp_star)) {
      perp_diff <- -Inf
    } else {
      perp_diff <- perp_star - perp_fixed
    }

    # Increasing one unit the number of tries
    tries <- tries + 1
  }
  return(list(rho = rho, prob_i = prob_star))
}

#' @title Binary search rho optimization for each observation
#'
#' @description Calculate the rho list values based on a fixed perplexity and a given data in \eqn{(\mathcal{S}^p)^r}.
#' The boundaries are set from 0 to 0.9999 for each value. It prints the time consumption.
#'
#' @inheritParams high_dimension
#' @inheritParams rho_optim_i_bst
#' @param cl defines a cluster to work with concurrently
#' @return Rho values and conditional probability matrix obtained as a result of the optimization.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(20, 3, 4)
#' rho_optim_bst(x, perp_fixed = 15, cl = clusterFactory(2))
#' rho_optim_bst(x, perp_fixed = 26, cl = clusterFactory(2))
rho_optim_bst <- function(x, perp_fixed,
                          cl = clusterFactory(parallel::detectCores() - 1)) {
  # Sample size
  n <- nrow(x)
  # Time start
  start_time <- Sys.time()
  # Concurrently calculate the optimized value of each observation
  res_opt <- parallel::parLapply(cl, 1:n, function(i) {
    rho_optim_i_bst(x, i, perp_fixed)
  })
  # Time end
  end_time <- Sys.time()
  # Print time difference
  print(end_time - start_time)
  # Stop clusters
  parallel::stopCluster(cl)
  # Aggregate the rho values of each observation in a list
  rho_values <- sapply(1:length(res_opt), function(i) res_opt[[i]]$rho)
  # Aggregate all the conditional probabilities given each observation in a matrix
  P <- t(sapply(1:length(res_opt), function(i) res_opt[[i]]$prob_i))
  return(list(rho_values = rho_values, P = P))
}

# > rho_optim_bst(x_array, 22)
# Time difference of 45.21138 secs
# with optim it took 2.58 min in the best case


# > rho_optim_bst(spokes, 22)
# Time difference of 5.758403 mins
# with optim func it took 28 min approx.
