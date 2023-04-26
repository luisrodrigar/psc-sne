###############################################
##          Polyspherical Cauchy HD          ##
## High-dimension neighborhood probabilities ##
###############################################

#' @title Polyspherical Cauchy conditional probability matrix (matrix version)
#'
#' @description Calculates the high-dimension conditional probabilities of a polyspherical Cauchy distribution. Matrix version algorithm.
#'
#' @param x an array of size \code{c(n, d + 1, r)} with the polyspherical data, where \code{n} is the number of observations, \code{d} is the dimension of each sphere, and \code{r} is the number of spheres.
#' @param rho_list rho list of size \code{n} for each \code{i}-th observation that stands for the concentration parameter.
#' @param cos_psh cosine similarities array of dimension \code{c(n, n, r)} for the polysphere.
#' @return An array of size \code{c(n, n)} with the high-dimension conditional probabilities of \code{x}.
#' @export
#' @examples
#' n <- 100
#' rho_list <- rep(seq(1, 10), each = n / 10)
#' x <- sphunif::r_unif_sph(n, 3, 3)
#' high_dimension_mat(x, rho_list)
#' high_dimension_mat(x, rho_list, cosine_polysph(x))
high_dim_mat <- function(x, rho_list, cos_psh = NULL) {

  if (!rlang::is_vector(rho_list)) {
    stop("rho_list must be a vector")
  }
  if (length(rho_list) != nrow(x)) {
    stop("rho_list size has to be equal to nrow(x)")
  }
  if (!is.null(cos_psh) && length(dim(cos_psh)) != 3) {
    stop("cos_sim_psh must be an array of size c(n, p + 1, r), from (S^p)^r")
  }
  # Number of observations
  n <- nrow(x)
  # Number of columns
  p <- ncol(x) - 1
  # Number of spheres
  r <- dim(x)[3]

  # Calculate the cosine similarities of 'x' if 'cos_sim_psh' param is null
  if (is.null(cos_psh)) {
    cos_psh <- sphunif::Psi_mat(x, scalar_prod = TRUE)
  }


  cos_psh_3d <- sapply(X = seq_len(r),
                       FUN = function(k) vec2matrix(cos_psh[, k], n = n, diag_value = 0),
                       simplify = 'array')

  # Calculate (1 + -2 * rho_list * (Y[i,,] %*% Y[j,,]) + (rho_list^2))^(-p) by each row of the 3d-array
  P <- (1 + rho_list^2 + cos_psh_3d * -2 * rho_list)^(-p)

  # Set the diagonal value to 0
  for (k in seq_len(r)) {

    P[, , k] = diag_3d(P, k, 0)

  }

  # Product operator by matrices of the 3d-array
  P_i_r <- apply(P, MARGIN = c(1,2), prod)
  # Summation operator by rows
  Pi <- rowSums(P_i_r)
  # Calculate (P_ij)_{ij} / (P_i)_{i}
  P_ij <- sweep(
    x = P_i_r, MARGIN = c(1, 2), STATS = Pi, FUN = "/")

  return(P_ij)
}


#' @title Polyspherical Cauchy conditional probability matrix
#'
#' @description Calculates the high-dimension conditional probabilities of a
#' polyspherical Cauchy distribution.
#'
#' @param x an array of size \code{c(n, d + 1, r)} with the polyspherical data,
#' where \code{n} is the number of observations, \code{d} is the dimension of
#' each sphere, and \code{r} is the number of spheres.
#' @param rho_list rho list of size \code{n} for each \eqn{i}-th observation
#' that stands for the concentration parameter.
#' @param cos_sim_psh a vector of size \code{n/2} with the cosine
#' similarities in high-dimension for the polysphere \eqn{(\mathcal{S}^p)^r}.
#' The way that the cosine similarities matrix is treated makes the calculus
#' faster since it is flat in a vector object. Optional parameter, defaults to
#' \code{NULL}.
#' @return An array of size \code{c(n, n)} with the high-dimension conditional
#' probabilities of \code{x}.
#' @inheritParams rho_optim_bst
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(100, 3, 3)
#' high_dimension(x, rep(0.5, 100), num_cores = 2)
#' high_dimension(x, rep(0.5, 100), sphunif::Psi_mat(x, scalar_prod = TRUE),
#'                num_cores = 2)
high_dimension <- function(x, rho_list, cos_sim_psh = NULL,
                           num_cores = parallel::detectCores() - 1) {
  if (!rlang::is_vector(rho_list)) {
    stop("rho_list must be a vector")
  }
  if (length(rho_list) != nrow(x)) {
    stop("rho_list size has to be equal to nrow(x)")
  }
  # Sample size
  n <- nrow(x)

  # Calculate the cosine similarities of 'x' if 'cos_sim_psh' param is null
  if (is.null(cos_sim_psh)) {
    cos_sim_psh <- sphunif::Psi_mat(x, scalar_prod = TRUE)
  }

  Pis <- parallel::parLapply(clusterFactory(num_cores), seq_len(n),
                             function(i) {
                               high_dimension_i(x = x, i = i, rho = rho_list[i],
                                                cos_sim_psh = cos_sim_psh)
  })

  return(do.call(rbind, Pis))
}


#' @title Polyspherical \eqn{i}-th Cauchy conditional probability vector
#' (matrix version)
#'
#' @description Calculates the high-dimension conditional probabilities of a
#' polyspherical Cauchy distribution for the i-th observation.
#' Matrix version algorithm.
#'
#' @inheritParams high_dimension
#' @inheritParams d_sph_cauchy
#' @return An vector of size \code{n} with the high-dimension conditional
#' probabilities of the \eqn{i}-th observation in \code{x}.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(100, 3, 3)
#' high_dimension_i(x, 1, 0.5)
#' high_dimension_i(x, 100, 0.5, sphunif::Psi_mat(x, scalar_prod = TRUE))
high_dimension_i <- function(x, i, rho, cos_sim_psh = NULL) {

  # sample size
  n <- nrow(x)
  # Dimension of each sphere
  p <- ncol(x) - 1
  # number of spheres
  r <- dim(x)[3]

  # Calculate the cosine similarities of 'x' if 'cos_sim_psh' param is null
  if (is.null(cos_sim_psh)) {
    cos_sim_psh <- sphunif::Psi_mat(x, scalar_prod = TRUE)
  }

  # Create vector of each cosine similarity for the i-th observation
  # with respect all the other observations
  cos_sim_i <- cos_sim_i(cos_sim_psh, i, r, n)

  # Applying formula of the polyspherical Cauchy density
  P <- (1 + rho^2 - 2 * rho * cos_sim_i)^(-p)
  # Set zero to the i-th observation
  P[i, ] <- 0

  # Product operator by matrices of the 3d-array
  P_i_r <- apply(P, MARGIN = 1, prod)

  # Summation operator by rows
  Pi <- sum(P_i_r)
  # Calculate (P_ij)_{ij} / (P_i)_{i}
  return(P_i_r / Pi)

}

