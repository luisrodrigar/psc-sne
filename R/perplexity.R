###########################################
##      spherical Cauchy Perplexity      ##
###########################################

#' @title Perplexity for the i-th observation (vector version)
#'
#' @description Calculate the perplexity of the i-th observation for a given
#' \eqn{\mathbf{\rho}} parameters list. Vector version algorithm.
#'
#' @inheritParams high_dimension
#' @inheritParams d_sph_cauchy
#' @return Perplexity of the \eqn{i}-th observation for all
#' the remainder observations.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(25, 3, 3)
#' i <- 1
#' to_perp_i(x, i, 0.25)
to_perp_i <- function(x, i, rho, cos_sim_psh = NULL) {
  # Calculate the cosine similarities of 'x' if 'cos_sim_psh' param is null
  if (is.null(cos_sim_psh)) {
    cos_sim_psh <- sphunif::Psi_mat(x, scalar_prod = TRUE)
  }
  # Obtaining the high dimension probability vector for the i-th observation
  P_i <- high_dimension_i(x, i, rho, cos_sim_psh)
  # Apply the formula of the entropy (matrix way)
  op <- P_i * log2(P_i)
  # Return the perplexity list result after applying the
  # perplexity formula (matrix way)
  return(list(perp = 2^(-sum(op, na.rm = TRUE)), Pjcondi = P_i))
}

#' @title Perplexity matrix (matrix version)
#'
#' @description Calculate the perplexity of each observations for a given
#' \eqn{\mathbf{\rho}} parameters list. Matrix version algorithm.
#'
#' @inheritParams high_dimension
#' @return Perplexity of each \eqn{i}-th observation for all the remainder
#' observations.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(25, 3, 3)
#' to_perp(x, rep(0.5, 25))
#' to_perp(x, rep(1 - 1e-4, 25), sphunif::Psi_mat(x, scalar_prod = TRUE))
to_perp <- function(x, rho_list, cos_sim_psh = NULL) {
  if (length(dim(x)) != 3) {
    stop("x must be an array of dimension c(n, p + 1, r), from (S^p)^r")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("rho_list must be a vector")
  }
  # Calculate the cosine similarities of 'x' if 'cos_sim_psh' param is null
  if (is.null(cos_sim_psh)) {
    cos_sim_psh <- sphunif::Psi_mat(x, scalar_prod = TRUE)
  }
  # Calculate the high-dimension probabilities
  P <- high_dimension(x, rho_list, cos_sim_psh)
  # Apply the formula of the entropy (matrix way)
  op <- P * log2(P)
  # Set the diagonal elements to zero
  diag(op) <- 0
  # Return the perplexity list result after applying the perplexity
  # formula (matrix way)
  return(2^(-rowSums(op)))
}


#####################################################
##  Optimize rho list based on a fixed perplexity  ##
#####################################################


#' @title Binary search rho optimization for each observation
#'
#' @description Calculate the rho list values based on a fixed perplexity and
#' a given data in \eqn{(\mathcal{S}^p)^r}. The boundaries are set from 0 to
#' 0.9999 for each value. It prints the time consumption.
#'
#' @inheritParams high_dimension
#' @inheritParams rho_optim_i_bst
#' @param num_cores number of cores to execute the code concurrently. This
#' value must be below the total number of the CPU has available.
#' @return Rho values and conditional probability matrix obtained as a
#' result of the optimization.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(20, 3, 4)
#' rho_optim_bst(x, perp_fixed = 15, num_cores = 2)
#' rho_optim_bst(x, perp_fixed = 26, num_cores = 2,
#'               sphunif::Psi_mat(x, scalar_prod = TRUE))
rho_optim_bst <- function(x, perp_fixed,
                          num_cores = parallel::detectCores() - 1,
                          cos_sim_psh = NULL) {
  if (is.null(cos_sim_psh)) {

    # Calculate the cosine similarities of (S^p)^r
    cos_sim_psh <- sphunif::Psi_mat(x, scalar_prod = TRUE)

  }
  # Sample size
  n <- nrow(x)
  # Time start
  start_time <- Sys.time()
  # Setting up the characteristics of the parallelization
  cl <- clusterFactory(num_cores)
  # Concurrently calculate the optimized value of each observation
  res_opt <- parallel::parLapply(cl, seq_len(n), function(i) {
    rho_optim_i_bst(x = x, i = i, perp_fixed = perp_fixed,
                    cos_sim_psh = cos_sim_psh)
  })
  # Time end
  end_time <- Sys.time()
  # Print time difference
  print(end_time - start_time)
  # Stop clusters
  parallel::stopCluster(cl)
  # Aggregate the rho values of each observation in a list
  rho_values <- sapply(seq_along(res_opt), function(i) {
    res_opt[[i]]$rho
  })
  # Aggregate all the conditional probabilities
  # given each observation in a matrix
  P <- t(sapply(seq_along(res_opt), function(i) {
    res_opt[[i]]$prob_i
  }))
  return(list(rho_values = rho_values, P = P))
}


#' @title Binary search rho optimization for the \eqn{i}-th observation
#'
#' @description  Calculate the rho based on a fixed perplexity
#' and a given data in \eqn{(\mathcal{S}^p)^r}. The boundaries are set
#' from 0 to 0.9999.
#'
#' @inheritParams high_dimension
#' @param i the \eqn{i}-th observation.
#' @param perp_fixed a fixed value used to optimized the rho values.
#' @param tolerance whether the difference between previous and current
#' results is below this value (optional, default value \code{1e-3}).
#' @param rho parameter which determines the concentration of the spherical
#' Cauchy distribution (optional, default value \code{0.5}).
#' @param max_tries number of maximum tries for each value of the rho list.
#' @return Rho concentration parameter and conditional probabilities calculated
#' for the \eqn{i}-th observation.
#' @keywords internal
rho_optim_i_bst <- function(x, i, perp_fixed, cos_sim_psh, tolerance = 1e-3,
                            rho = 0.5, max_tries = 20) {
  # Min value of rho
  rho_min <- 0
  # Max value of rho
  rho_max <- 1 - 1e-4
  # Number of tries
  tries <- 0
  # Difference of perplexity, initially set to NA
  perp_diff <- NA

  # Calculate perplexity based on the default value of rho
  res <- to_perp_i(x = x, i = i, rho = rho, cos_sim_psh = cos_sim_psh)
  # Get the perplexity value
  perp_star <- res$perp
  # Get the conditional probabilities given i
  prob_star <- res$Pjcondi

  # Set the difference to - infinity if perplexity is NA,
  # calculate difference otherwise
  if (is.na(perp_star)) {
    perp_diff <- -Inf
  } else {
    perp_diff <- perp_star - perp_fixed
  }

  # While difference of perplexity is NA or greater than tolerance
  # and the number of tries less than the maximum, do the following
  while ((is.na(perp_diff) || abs(perp_diff) > tolerance) &&
         tries <= max_tries) {

    # Find by means of the binary search optimal rho
    if (perp_diff > 0 || is.na(perp_diff)) {
      rho_min <- rho
    } else {
      rho_max <- rho
    }
    rho <- (rho_min + rho_max) / 2

    # Calculate the perplexity and the P_j|i for the values of i-th iteration
    res_loop <- to_perp_i(x = x, i = i, rho = rho, cos_sim_psh = cos_sim_psh)
    perp_star <- res_loop$perp
    prob_star <- res_loop$Pjcondi

    # Again, calculate the difference of perplexity as it was done
    # outside the loop
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

