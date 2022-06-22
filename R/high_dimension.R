###############################################
##        Poly-spherical Cauchy HD           ##
## High-dimension neighborhood probabilities ##
###############################################

library(lsa)


#' Calculate the symmetric probabilities of a given conditional probabilities matrix
#'
#' @param P matrix of probabilities (P_i|j)_{ij}
#' @return The sum of \code{P} and \code{t(P)} divided by twice the sample size
#' @examples
#' symmetric_probs(matrix(runif(3 * 3), nrow = 3, ncol = 3))
#' symmetric_probs(diag(3))
symmetric_probs <- function(P) {
  n <- nrow(P)
  P <- (P + t(P)) / (2 * n)
  return(P)
}

#' Set a specific value to the diagonal of a matrix
#'
#' @param x array 3-dimensional
#' @param k the k-th matrix of the 3d array
#' @param val the specific value to set in the diagonal
#' @return The \code{k}-th matrix of \code{x} with the diagonal set to \code{val}
#' @examples
#' diag_3d(x, 1, 0)
#' diag_3d(x, 3, 0)
diag_3d <- function(x, k, val) {
  if (k < 1 || k > dim(x)[3]) {
    stop("The 3rd dimensional index k not valid, must be >= 1 and <= r")
  }
  diag(x[, , k]) <- val
  return(x[, , k])
}

#' Calculates the cosine similarity for each matrix
#'
#' @param x array 3-dimensional
#' @return a 3-dimensional array with the cosine similarities of each matrix of \code{x}
#' @examples
#' cosine_polysph(x)
cosine_polysph <- function(x) {
  r <- dim(x)[3]
  cosine_sphere_ith <- function(k) {
    cos_sim <- cosine(t(x[, , k]))
  }
  sapply(1:r, cosine_sphere_ith, simplify = "array")
}

#' Matrix version
#' Calculates the high-dimension probabilities of a polyspherical cauchy distribution
#'
#' @param x array 3-dimensional
#' @param rho_list rho parameter for each \code{i}-th observation
#' @param cos_sim_pol cosine similarities of the polysphere
#' @return a 3-dimensional array with the high-dimension probabilities of the array \code{x}
#' @examples
#' high_dimension(x, rep(0.5, 3))
#' high_dimension(x, rep(0.5, 3), cosine_polysph(diag(3)))
high_dimension <- function(x, rho_list, cos_sim_pol = NULL) {
  if (!rlang::is_vector(rho_list)) {
    stop("Parameter rho_list must be a vector")
  }
  if (length(rho_list) != nrow(x)) {
    stop("Parameter rho_list not valid, size must be equal to nrow(x)")
  }
  if (!is.null(cos_sim_pol) && length(dim(cos_sim_pol)) != 3) {
    stop("Parameter cos_sim_ps must be a 3d-array")
  }
  # Number of observations
  n <- nrow(x)
  # Number of columns
  p <- ncol(x) - 1
  # Number of spheres
  r <- dim(x)[3]

  # Calculate the cosine similarities of 'x' if 'cos_sim_pol' param is null
  if (is.null(cos_sim_pol)) {
    cos_sim_pol <- cosine_polysph(x)
  }

  # Calculate -2 * rho_list * (Y[i,,] %*% Y[j,,]) by each row of the 3d-array
  P <- sweep(cos_sim_pol,
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

  # Set the diagonal of each matrix of the 3d-array to zero
  Paux <- P
  P <- sapply(1:r, FUN = diag_3d, x = Paux, val = 0, simplify = "array")
  # Product operator by matrices of the 3d-array
  P_i_r <- apply(P, MARGIN = c(1, 2), prod)
  # Summation operator by rows
  Pi <- rowSums(P_i_r)
  # Calculate (P_ij)_{ij} / (P_i)_{i}
  P_ij <- sweep(
    x = P_i_r, MARGIN = c(1, 2), STATS = Pi, FUN = "/",
    check.margin = FALSE
  )
  return(P_ij)
}

#' Scalar version
#' Calculates the high-dimension probabilities of a polyspherical Cauchy distribution
#'
#' @param x array 3-dimensional
#' @param rho_list rho parameter for each \code{i}-th observation
#' @return a 3-dimensional matrix with the high-dimension probabilities of the array \code{x}
#' @examples
#' high_dimension_p(x, rep(0.5, 3))
#' high_dimension_p(x, rep(0.5, 3))
high_dimension_p <- function(x, rho_list) {
  if (!rlang::is_vector(rho_list)) {
    stop("Parameter rho_list must be a vector")
  }
  if (length(rho_list) != nrow(x)) {
    stop("Parameter rho_list not valid, size must be equal to nrow(x)")
  }
  # Marginal conditional probability of the i-th observation
  total_p <- P_total_psc(x, rho_list)
  # Function that calculates the probability of the j-th observation given the i-th one
  jcondi <- function(i) {
    sapply(1:nrow(x), function(j) {
      jcondi_psc(x, i, j, rho_list, total_p[i])
    })
  }
  # Apply the previous function for each observation
  return(t(sapply(1:nrow(x), jcondi)))
}
