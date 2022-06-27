##############################################
##        Poly-spherical Cauchy HD          ##
## Low-dimension neighborhood probabilities ##
##############################################

#' Projection of the points onto the sphere of radius 1
#'
#' @param y matrix with the points in the sphere
#' @return a matrix with the values of x projected onto the sphere of radius 1
#' @examples
#' y <- rotasym::r_unif_sphere(100, 2)
#' radial_projection(y)
radial_projection <- function(y) {
  # For every observation, apply the formula of the radial projection
  t(sapply(1:nrow(y), function(i) {
    # y_i / |y_i|
    y[i, ] / norm(y[i, ], type = "2")
  }))
}

#' Calculate the low-dimension probabilities of a reduced matrix Y
#' Using
#'
#' @param Y matrix with the points in the sphere
#' @param rho parameter between 0 and 1 (not included)
#' @return a matrix with the values of x projected onto the sphere of radius 1
#' @examples
#' Y <- rotasym::r_unif_sphere(100, 2)
#' low_dimension_Q(Y, 0)
#' low_dimension_Q(Y, 0.5)
#' low_dimension_Q(Y, 0.9999)
low_dimension_Q <- function(Y, rho) {
  # Obtaining d, where S^d
  d <- ncol(Y) - 1
  # Projecting the points onto de sphere, in case they are not
  Z <- radial_projection(Y)
  # Calculate the cosine similarities matrix of Z
  cos_simil <- lsa::cosine(t(Z))
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
