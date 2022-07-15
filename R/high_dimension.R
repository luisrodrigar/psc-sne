###############################################
##        Polyspherical Cauchy HD           ##
## High-dimension neighborhood probabilities ##
###############################################

#' @title Polyspherical Cauchy conditional probability matrix (matrix version)
#'
#' @description Calculates the high-dimension conditional probabilities of a polyspherical Cauchy distribution. Matrix version algorithm.
#'
#' @param x an array of size \code{c(n, d + 1, r)} with the polyspherical data, where \code{n} is the number of observations, \code{d} is the dimension of each sphere, and \code{r} is the number of spheres.
#' @param rho_list rho list of size \code{n} for each \code{i}-th observation that stands for the concentration parameter.
#' @param cos_sim_psh cosine similarities array of dimension \code{c(n, n, r)} for the polysphere.
#' @return An array of size \code{c(n, n)} with the high-dimension conditional probabilities of \code{x}.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(100, 3, 3)
#' high_dimension(x, rep(0.5, 100))
#' high_dimension(x, rep(0.5, 100), cosine_polysph(x))
high_dimension <- function(x, rho_list, cos_sim_psh = NULL) {
  if (!rlang::is_vector(rho_list)) {
    stop("rho_list must be a vector")
  }
  if (length(rho_list) != nrow(x)) {
    stop("rho_list size has to be equal to nrow(x)")
  }
  if (!is.null(cos_sim_psh) && length(dim(cos_sim_psh)) != 3) {
    stop("cos_sim_psh must be an array of size c(n, p + 1, r), from (S^p)^r")
  }
  # Number of observations
  n <- nrow(x)
  # Number of columns
  p <- ncol(x) - 1
  # Number of spheres
  r <- dim(x)[3]

  # Calculate the cosine similarities of 'x' if 'cos_sim_psh' param is null
  if (is.null(cos_sim_psh)) {
    cos_sim_psh <- sphunif::Psi_mat(x, scalar_prod = TRUE)
  }

  # Calculate -2 * rho_list * (Y[i,,] %*% Y[j,,]) by each row of the 3d-array
  P <- sweep(cos_sim_psh,
             MARGIN = 1, STATS = (-2 * rho_list), FUN = "*",
             check.margin = FALSE
  )
  # Calculate P + (rho_list^2) by each row of the 3d-array
  P <- sweep(P,
             MARGIN = 1, STATS = (rho_list^2), FUN = "+",
             check.margin = FALSE
  )
  # Calculate 1 / (1 + P)^p
  P <- 1 / (1 + P)^p

  # Product operator by matrices of the 3d-array
  P_i_r <- apply(P, MARGIN = 1, prod)
  # Reconstruct from vector to symmetric matrix
  P_i_r <- vec2matrix(P_i_r, n, diag_value = 0)

  # Summation operator by rows
  Pi <- rowSums(P_i_r)
  # Calculate (P_ij)_{ij} / (P_i)_{i}
  P_ij <- sweep(
    x = P_i_r, MARGIN = c(1, 2), STATS = Pi, FUN = "/",
    check.margin = FALSE
  )
  return(P_ij)
}

#' @title Polyspherical Cauchy conditional probability matrix (scalar version)
#'
#' @description Calculates the high-dimension conditional probabilities of a polyspherical Cauchy distribution. Scalar version algorithm.
#'
#' @inheritParams high_dimension
#' @return An array of size \code{c(n, n)} with the high-dimension conditional probabilities of \code{x}.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(100, 3, 3)
#' high_dimension_p(x, rep(0.5, 100))
high_dimension_p <- function(x, rho_list) {
  if (!rlang::is_vector(rho_list)) {
    stop("rho_list must be a vector")
  }
  if (length(rho_list) != nrow(x)) {
    stop("rho_list size has to be equal to nrow(x)")
  }
  # Marginal conditional probability of the i-th observation
  d_total_psph <- d_total_psph_cauchy(x, rho_list)
  # Function that calculates the probability of the j-th observation given the i-th one
  jcondi <- function(i) {
    sapply(1:nrow(x), function(j) jcondi_psph(x, i, j, rho_list, d_total_psph[i]))
  }
  # Apply the previous function for each observation
  return(t(sapply(1:nrow(x), jcondi)))
}
