##############################################
##           Spherical Cauchy HD            ##
## Low-dimension neighborhood probabilities ##
##############################################

#' @title Low-dimension probabilities
#'
#' @description Calculate the low-dimension probabilities of a reduced matrix Y.
#'
#' @param Y matrix of size \code{c(n, d)}, where \code{n} is the number of observation, with the points onto the sphere \eqn{\mathcal{S}^d}.
#' @inheritParams d_sph_cauchy
#' @return A matrix with the values of x projected onto the sphere of radius 1.
#' @export
#' @examples
#' Y <- rotasym::r_unif_sphere(100, 2)
#' low_dimension_Q(Y, 0)
#' low_dimension_Q(Y, 0.5)
#' low_dimension_Q(Y, 0.9999)
low_dimension_Q <- function(Y, rho) {
  # Obtaining d, where S^d
  d <- ncol(Y) - 1
  # Obtaining the sample size
  n <- nrow(Y)
  # Projecting the points onto the sphere, in case they are not
  Z <- radial_projection(Y)
  # Calculate the cosine similarities matrix of Z
  cos_simil <- sphunif::Psi_mat(array(Z, dim = c(n, d + 1, 1)), scalar_prod = TRUE)

  # Applying the Spherical Cauchy low-dimension joint function
  Q <- (1 + rho^2 - 2 * rho * cos_simil)^(-d)

  # Calculate the total probability
  Qi <- 2 * sum(Q)
  # Calculate the matrix probability of the joint probabilities
  Q_ij <- Q / Qi

  # Reconstruct from vector to symmetric matrix
  Q_ij <- vec2matrix(Q_ij, n, 0)

  return(Q_ij)
}


