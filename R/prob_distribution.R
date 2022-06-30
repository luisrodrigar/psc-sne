#####################################
### Spherical Cauchy Distribution ###
#####################################

####################
## High-dimension ##
####################

#' @title High-dimension polyspherical Cauchy density
#'
#' @description Calculate the high-dimension spherical Cauchy density function.
#'
#' @inheritParams high_dimension
#' @param i corresponds to the \eqn{i}-th observation index.
#' @param j corresponds to the \eqn{j}-th observation index.
#' @param rho concentration parameter must be in [0, 1).
#' @param k corresponds to the k-th sphere index.
#' @param p corresponds to S^p, associated to R^(p+1).
#' @return Spherical Cauchy density value given the parameters.
d_sph_cauchy <- function(x, i, j, rho, k, p) {
  if (i < 1 || i > nrow(x) || j < 1 || j > nrow(x)) {
    stop("i or j not valid, must be within 1 and nrow(x)")
  }
  if (k < 1 || k > dim(x)[3]) {
    stop("k index not valid, must be between 1 and r, where (S^p)^r")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("rho must be an scalar")
  }
  if (p + 1 != ncol(x)) {
    stop("p is not valid, it must match with the dataset dimension ncol(x)-1")
  }
  # Applying the formula: 1 + rho^2 - 2 * rho * x[i,,k]'x[j,,k]
  drop((1 + rho^2 - 2 * rho * t(x[i, , k]) %*% x[j, , k])^(-p))
}

#' @title Polyspherical Cauchy density for the \eqn{i}-th and \eqn{j}-th observations
#'
#' @description Calculate the high-dimension polyspherical Cauchy density for the
#' \eqn{i}-th and \eqn{j}-th observations.
#'
#' @inheritParams high_dimension
#' @inheritParams d_sph_cauchy
#' @return Polyspherical Cauchy density value given the parameters.
d_ij_psph_cauchy <- function(x, i, j, rho) {
  if (i < 1 || i > nrow(x) || j < 1 || j > nrow(x)) {
    stop("i or j not valid, must be between 1 and nrow(x)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("rho must be an scalar")
  }
  # If the i-th and the j-th observation are the same, the value returned is zero
  if (i == j) {
    return(0)
  }
  # Productory of the density of the i-th and j-th observations for each sphere (S^p)
  return(prod(sapply(1:(dim(x)[3]),
    FUN = d_sph_cauchy,
    x = x,
    i = i,
    j = j,
    rho = rho,
    p = (ncol(x) - 1)
  )))
}

#' @title Marginal polyspherical Cauchy density for the \eqn{i}-th observation

#' @description Calculate the high-dimension polyspherical Cauchy density for the
#' \eqn{i}-th observation.
#'
#' @inheritParams high_dimension
#' @inheritParams d_sph_cauchy
#' @return Polyspherical Cauchy density value given the parameters.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(20, 3, 4)
#' d_i_psph_cauchy(x, 1, rep(0.5, 20))
#' d_i_psph_cauchy(x, 4, rep(0.9999, 20))
d_i_psph_cauchy <- function(x, i, rho_list) {
  if (i < 1 || i > nrow(x)) {
    stop("i not valid, must be between 1 and nrow(x)")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("rho must be a vector of dimension n")
  }
  # Sum of the probabilities of the productory for a given i-th observation
  return(sapply(1:nrow(x), FUN = d_ij_psph_cauchy, x = x, i = i, rho = rho_list[i]))
}

#' @title Polyspherical marginal density function values for all the observations
#'
#' @description Calculate the marginal high-dimension polyspherical Cauchy probabilities for the \eqn{i}-th observation.
#'
#' @inheritParams high_dimension
#' @return Marginal polyspherical Cauchy probabilities vector given the \eqn{\boldsymbold{\rho}} parameters.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(20, 3, 4)
#' d_total_psph_cauchy(x, rep(0.5, 20))
#' d_total_psph_cauchy(x, rep(0.9999, 20))
d_total_psph_cauchy <- function(x, rho_list) {
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("rho must be a vector of dimension n")
  }
  # Sum of the row probabilities to get the marginal ones of each i-th observation
  return(rowSums(sapply(1:nrow(x),
    FUN = d_i_psph_cauchy, x = x, rho = rho_list,
    simplify = "array"
  )))
}

#' @title Conditional polyspherical Cauchy probability for the \eqn{i}-th and \eqn{j}-th observations
#'
#' @description Calculate the conditional high-dimension probability of the \eqn{j}-th observation given the \eqn{i}-th.
#'
#' @inheritParams high_dimension
#' @inheritParams d_sph_cauchy
#' @param d_total_i_psph_cauchy marginal probability of the \eqn{i}-th observation. Optional, defaulta value \code{NULL}.
#' @return Conditional polyspherical Cauchy probability of the \eqn{j}-th given the \eqn{i}-th observation.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(20, 3, 4)
#' jcondi_psph(x, 1, 2, rep(0.5, 20), d_total_psph_cauchy(x, rep(0.5, 20))[1])
#' jcondi_psph(x, 4, 6, rep(0.9999, 20), d_total_psph_cauchy(x, rep(0.9999, 20))[4])
jcondi_psph <- function(x, i, j, rho_list, d_total_i_psph_cauchy = NULL) {
  if (i < 1 || i > nrow(x) || j < 1 || j > nrow(x)) {
    stop("i or j not valid, must be between 1 and nrow(x)")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("rho must be a vector of dimension n")
  }
  # Calculate the marginal density for the i-th observation
  if (is.null(d_total_i_psph_cauchy)) {
    d_total_i_psph_cauchy <- sum(d_i_psph_cauchy(x, i, rho_list))
  }
  if (!rlang::is_scalar_atomic(d_total_i_psph_cauchy)) {
    stop("d_total_i_psph_cauchy must be an scalar")
  }
  # If the i-th and the j-th observation are the same, the value returned is zero
  if (i == j) {
    return(0)
  }
  # Conditional probability of the j-th given the i-th observation
  return(d_ij_psph_cauchy(x, i, j, rho_list[i]) / d_total_i_psph_cauchy)
}

#' @title Conditional polyspherical Cauchy probability for the \eqn{i}-th observation
#'
#' @description Calculate the conditional high-dimension probability for all the \eqn{j}-th observation given the \eqn{i}-th.
#'
#' @inheritParams high_dimension
#' @inheritParams d_sph_cauchy
#' @inheritParams jcondi_psph
#' @return Conditional polyspherical Cauchy probability for all the \eqn{j}-th given the \eqn{i}-th observation.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(20, 3, 4)
#' rho_list_1 <- rep(0.5, 20)
#' rho_list_2 <- rep(0.9999, 20)
#' d_total_1_psph <- d_total_psph_cauchy(x, rho_list_1)[1]
#' d_total_4_psph <- d_total_psph_cauchy(x, rho_list_2)[4]
#' prob_cond_i_psph(x, 1, rho_list_1, d_total_1_psph)
#' prob_cond_i_psph(x, 4, rho_list_2, d_total_4_psph)
prob_cond_i_psph <- function(x, i, rho_list, d_total_i_psph_cauchy = NULL) {
  if (i < 1 || i > nrow(x)) {
    stop("i not valid, must be between 1 and nrow(x)")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("rho_list must be a vector of dimension n")
  }
  # Calculate the marginal density for the i-th observation
  if (is.null(d_total_i_psph_cauchy)) {
    d_total_i_psph_cauchy <- sum(d_i_psph_cauchy(x, i, rho_list))
  }
  if (!rlang::is_scalar_atomic(d_total_i_psph_cauchy)) {
    stop("d_total_i_psph_cauchy must be an scalar")
  }
  # Calculate every j-th conditional probability given the i-th observation
  return(sapply(1:nrow(x), jcondi_psph,
    x = x, i = i, rho_list = rho_list,
    d_total_i_psph_cauchy = d_total_i_psph_cauchy
  ))
}

#####################
##  Low-dimension  ##
#####################

#' @title Marginal spherical Cauchy density function
#'
#' @description Calculate the marginal low-dimension spherical Cauchy probability.
#'
#' @param y matrix that stands for an sphere, \eqn{\mathcal{S}^d}.
#' @inheritParams d_sph_cauchy
#' @return Spherical marginal Cauchy density probability given the parameters.
d_i_sph_cauchy <- function(y, i, rho) {
  if (i < 1 || i > nrow(y)) {
    stop("i not valid, must be between 1 and nrow(y)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("rho must be an scalar")
  }
  # Sum of all the densities for the i-th observation and every all observations
  sum(sapply((1:nrow(y))[-i], d_sph_cauchy,
    x = array(y, dim = c(nrow(y), ncol(y), 1)),
    i = i, rho = rho, k = 1, p = ncol(y) - 1
  ))
}

#' @title Conditional spherical Cauchy probability
#'
#' @description Calculate the probability of choosing a pair of elements
#' where the \eqn{i}-th and the \eqn{j}-th observations are selected.
#'
#' @inheritParams d_i_sph_cauchy
#' @inheritParams d_sph_cauchy
#' @param d_i_sph_cauchy the marginal probability of the \eqn{i}-th observation. Optional, default value to \code{NULL}.
#' @return Spherical marginal Cauchy density probability given the parameters.
jcondi_sph <- function(y, i, j, rho, d_i_sph_cauchy = NULL) {
  if (i < 1 || i > nrow(y)) {
    stop("i not valid, must be between 1 and nrow(y)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("rho must be an scalar")
  }
  if (is.null(d_i_sph_cauchy)) {
    d_i_sph_cauchy <- d_i_sph_cauchy(y, i, rho)
  }
  # If the i-th and the j-th observation are the same, the value returned is zero
  if (i == j) {
    return(0)
  }
  # Calculate the probability of choosing one pairs at random (i-th and j-th observation)
  return(d_sph_cauchy(
    array(y, dim = c(nrow(y), ncol(y), 1)), i, j, rho, 1,
    ncol(y) - 1
  ) / d_i_sph_cauchy)
}

#' @title Conditional spherical Cauchy probabilities given the \eqn{i}-th observation
#'
#' @description Calculate \eqn{n}, where it stands for the sample size, probabilities where each element is the probability of choosing the
#' \eqn{i}-th and the \eqn{j}-th observations, the \eqn{i}-th element is a fixed element in the chosen pairs.
#'
#' @inheritParams d_i_sph_cauchy
#' @inheritParams d_sph_cauchy
#' @return Spherical marginal Cauchy density probability given the parameters.
#' @export
#' @examples
#' y <- rotasym::r_unif_sphere(100, 2)
#' prob_cond_i_sph(y, 1, 0.5)
#' prob_cond_i_sph(y, 4, 0.9999)
prob_cond_i_sph <- function(y, i, rho) {
  if (i < 1 || i > nrow(y)) {
    stop("i not valid, must be between 1 and nrow(y)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("rho must be an scalar")
  }
  # Calculate the probabilities of chosen a pair of elements where the i-th element is fixed
  return(sapply(1:nrow(y), jcondi_sph, y = y, i = i, rho))
}
