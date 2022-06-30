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
  cos_simil <- reconstruct_cos_sim_mat(
    cos_sim_vec = sphunif::Psi_mat(array(Z, dim = c(n, d + 1, 1)), scalar_prod = TRUE),
    n = n
  )
  # Applying the Spherical Cauchy low-dimension joint function
  Q <- (1 + rho^2 - 2 * rho * cos_simil)^(-d)
  # Set the elements of the diagonal to zero
  diag(Q) <- 0
  # Calculate the total probability
  Qi <- sum(Q)
  # Calculate the matrix probability of the joint probabilities
  Q_ij <- Q / Qi
  return(Q_ij)
}
