#####################################
### Spherical Cauchy Distribution ###
#####################################

####################
## High-dimension ##
####################

#' Calculate the high-dimension spherical Cauchy density function
#'
#' @param x 3d-array that is a poly-sphere, (S^p)^r
#' @param i corresponds to the i-th observation index
#' @param j corresponds to the j-th observation index
#' @param rho param between 0 and 1 (not included)
#' @param k corresponds to the k-th sphere index
#' @param p corresponds to S^p, associated to R^(p+1)
#' @return spherical Cauchy density value given the parameters
#' @examples
#' x_1 <- gen_polysphere(20, 1, 4)
#' simple_dspcauchy_hd(x_1, 4, 3, 0.9999, 3, 1)
#' x_2 <- gen_polysphere(20, 2, 4)
#' simple_dspcauchy_hd(x_2, 1, 2, 0.5, 2, 2)
simple_dspcauchy_hd <- function(x, i, j, rho, k, p) {
  if (i < 1 || i > nrow(x) || j < 1 || j > nrow(x)) {
    stop("The indexes (i and j) not valid, must be >= 1 and <= nrow(x)")
  }
  if (k < 1 || k > dim(x)[3]) {
    stop("The 3rd dimensional index k not valid, must be >= 1 and <= r")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("Parameter rho must be an scalar")
  }
  if (p + 1 != ncol(x)) {
    stop("Parameter p is not valid, it must match with the dataset dimension")
  }
  # Applying the formula: 1 + rho^2 - 2 * rho * x[i,,k]'x[j,,k]
  drop((1 + rho^2 - 2 * rho * t(x[i, , k]) %*% x[j, , k])^(-p))
}

#' Calculate the high-dimension poly-spherical Cauchy density for the
#' i-th and j-th observations.
#'
#' @param x 3d-array that is a poly-sphere, (S^p)^r
#' @param i corresponds to the i-th observation index
#' @param j corresponds to the j-th observation index
#' @param rho param between 0 and 1 (not included)
#' @return poly-spherical Cauchy density value given the parameters
#' @examples
#' x <- gen_polysphere(20, 2, 4)
#' P_ij_psc(x, 1, 2, 0.5)
#' P_ij_psc(x, 4, 3, 0.9999)
P_ij_psc <- function(x, i, j, rho) {
  if (i < 1 || i > nrow(x) || j < 1 || j > nrow(x)) {
    stop("The indexes (i and j) not valid, must be >= 1 and <= nrow(x)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("Parameter rho must be an scalar")
  }
  # If the i-th and the j-th observation are the same, the value returned is zero
  if (i == j) {
    return(0)
  }
  # Productory of the density of the i-th and j-th observations for each sphere (S^p)
  return(prod(sapply(1:(dim(x)[3]),
    simple_dspcauchy_hd,
    x = x, i = i, j = j, rho = rho, p = (ncol(x) - 1)
  )))
}

#' Calculate the high-dimension poly-spherical Cauchy density for the
#' i-th observation.
#'
#' @param x 3d-array that is a poly-sphere, (S^p)^r
#' @param i corresponds to the i-th observation index
#' @param rho_list param values between 0 and 1 (not included)
#' @return poly-spherical Cauchy density value given the parameters
#' @examples
#' x <- gen_polysphere(20, 2, 4)
#' P_i_psc(x, 1, rep(0.5, 20))
#' P_i_psc(x, 4, rep(0.9999, 20))
P_i_psc <- function(x, i, rho_list) {
  if (i < 1 || i > nrow(x)) {
    stop("The indexes i not valid, must be >= 1 and <= nrow(x)")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("Parameter rho must be a vector of dimension n")
  }
  # Sum of the probabilities of the productory for a given i-th observation
  return(sapply(1:nrow(x), FUN = P_ij_psc, x = x, i = i, rho = rho_list[i]))
}

#' Calculate the marginal high-dimension poly-spherical Cauchy probabilities for the i-th observation.
#'
#' @param x 3d-array that is a poly-sphere, (S^p)^r
#' @param rho_list param values between 0 and 1 (not included)
#' @return marginal poly-spherical Cauchy probabilities value given the parameters
#' @examples
#' x <- gen_polysphere(20, 2, 4)
#' P_total_psc(x, rep(0.5, 20))
#' P_total_psc(x, rep(0.9999, 20))
P_total_psc <- function(x, rho_list) {
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("Parameter rho must be a vector of dimension n")
  }
  # Sum of the row probabilities to get the marginal ones of each i-th observation
  return(rowSums(sapply(1:nrow(x),
    FUN = P_i_psc, x = x, rho = rho_list,
    simplify = "array"
  )))
}

#' Calculate the conditional high-dimension probability of the j-th observation given the i-th.
#'
#' @param x 3d-array that is a poly-sphere, (S^p)^r
#' @param i corresponds to the i-th observation index
#' @param j corresponds to the j-th observation index
#' @param rho_list param values between 0 and 1 (not included)
#' @param total_P_i marginal probability of the i-th observation
#' @return conditional poly-spherical Cauchy probability of the j-th given the i-th observation
#' @examples
#' x <- gen_polysphere(20, 2, 4)
#' jcondi_psc(x, 1, 2, rep(0.5, 20), P_total_psc(x, rep(0.9999, 20))[1])
#' jcondi_psc(x, 4, 6, rep(0.9999, 20), P_total_psc(x, rep(0.9999, 20))[4])
jcondi_psc <- function(x, i, j, rho_list, total_P_i) {
  if (i < 1 || i > nrow(x) || j < 1 || j > nrow(x)) {
    stop("The indexes (i and j) not valid, must be >= 1 and <= nrow(x)")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("Parameter rho must be a vector of dimension n")
  }
  if (!rlang::is_scalar_atomic(total_P_i)) {
    stop("Parameter total_P_i must be an scalar")
  }
  # If the i-th and the j-th observation are the same, the value returned is zero
  if (i == j) {
    return(0)
  }
  # Conditional probability of the j-th given the i-th observation
  return(P_ij_psc(x, i, j, rho_list[i]) / total_P_i)
}

#' Calculate the conditional high-dimension probability for all the j-th observation given the i-th.
#'
#' @param x 3d-array that is a poly-sphere, (S^p)^r
#' @param i corresponds to the i-th observation index
#' @param rho_list param values between 0 and 1 (not included)
#' @param total_P_i marginal probability of the i-th observation
#' @return conditional poly-spherical Cauchy probability for all the j-th given the i-th observation
#' @examples
#' x <- gen_polysphere(20, 2, 4)
#' psc_cond_given_i(x, 1, rep(0.5, 20), P_total_psc(x, rep(0.5, 20))[1])
#' psc_cond_given_i(x, 4, rep(0.9999, 20), P_total_psc(x, rep(0.9999, 20))[4])
psc_cond_given_i <- function(x, i, rho_list, total_P_i = NULL) {
  if (i < 1 || i > nrow(x)) {
    stop("The index i not valid, must be >= 1 and <= nrow(x)")
  }
  if (rlang::is_scalar_atomic(rho_list)) {
    stop("Parameter rho must be a vector of dimension n")
  }
  if (!rlang::is_scalar_atomic(total_P_i)) {
    stop("Parameter total_P_i must be an scalar")
  }
  # Calculate the marginal probability if it is null the parameter
  if (is.null(total_P_i)) {
    total_P_i <- sum(P_i_psc(x, i, rho_list))
  }
  # Calculate every j-th conditional probability given the i-th observation
  return(sapply(1:nrow(x), jcondi_psc,
    x = x, i = i, rho_list = rho_list,
    total_P_i = total_P_i
  ))
}

#####################
##  Low-dimension  ##
#####################

#' Calculate the low-dimension spherical Cauchy density function
#'
#' @param y matrix that stands for an sphere, S^d
#' @param i corresponds to the i-th observation index
#' @param j corresponds to the j-th observation index
#' @param rho param between 0 and 1 (not included)
#' @param d corresponds to S^d, associated to R^(d+1)
#' @return spherical Cauchy density value given the parameters
#' @examples
#' y_1 <- rotasym::r_unif_sphere(100, 2)
#' simple_dspcauchy_ld(y_1, i = 1, j = 2, rho = 0.5, d = 1)
#' y_2 <- rotasym::r_unif_sphere(100, 3)
#' simple_dspcauchy_ld(y_2, i = 4, j = 3, rho = 0.9999, d = 2)
simple_dspcauchy_ld <- function(y, i, j, rho, d) {
  if (i < 1 || i > nrow(y) || j < 1 || j > nrow(y)) {
    stop("The indexes (i and j) not valid, must be >= 1 and <= nrow(y)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("Parameter rho must be an scalar")
  }
  if (d + 1 != ncol(y)) {
    stop("Parameter p is not valid, it must match with the dataset dimension")
  }
  drop((1 + rho^2 - 2 * rho * t(y[i, ]) %*% y[j, ])^(-d))
}

#' Calculate the marginal low-dimension spherical Cauchy probability
#'
#' @param y matrix that stands for an sphere, S^d
#' @param i corresponds to the i-th observation index
#' @param rho param between 0 and 1 (not included)
#' @return spherical marginal Cauchy density probability given the parameters
#' @examples
#' y <- rotasym::r_unif_sphere(100, 2)
#' Q_i_sc(y, 1, 0.5)
#' Q_i_sc(y, 4, 0.9999)
Q_i_sc <- function(y, i, rho) {
  if (i < 1 || i > nrow(y)) {
    stop("The index i not valid, must be >= 1 and <= nrow(y)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("Parameter rho must be an scalar")
  }
  # Sum of all the densities for the i-th observation and every all observations
  sum(sapply((1:nrow(y))[-i], simple_dspcauchy_ld, y = y, i = i, rho = rho, d = ncol(y) - 1))
}

#' Calculate the probability of choosing a pair of elements
#' where the i-th and the j-th observations are selected
#'
#' @param y matrix that stands for a sphere, S^d
#' @param i corresponds to the i-th observation index
#' @param j corresponds to the j-th observation index
#' @param rho param between 0 and 1 (not included)
#' @param total_Q_i is the marginal probability of the i-th observation
#' @return spherical marginal Cauchy density probability given the parameters
#' @examples
#' y <- rotasym::r_unif_sphere(100, 2)
#' jcondi_sc(y, 1, 3, 0.5)
#' jcondi_sc(y, 4, 5, 0.9999, Q_i_sc(y, 4, 0.9999))
jcondi_sc <- function(y, i, j, rho, total_Q_i = NULL) {
  if (i < 1 || i > nrow(y)) {
    stop("The index i not valid, must be >= 1 and <= nrow(y)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("Parameter rho must be an scalar")
  }
  if (is.null(total_Q_i)) {
    total_Q_i <- Q_i_sc(y, i, rho)
  }
  # If the i-th and the j-th observation are the same, the value returned is zero
  if (i == j) {
    return(0)
  }
  # Calculate the probability of choosing one pairs at random (i-th and j-th observation)
  return(simple_dspcauchy_ld(y, i, j, rho, ncol(y) - 1) / total_Q_i)
}

#' Calculate n probabilities where each element is the probability of choosing the
#' i-th and the j-th observations, the i-th element is a fixed element in the chosen pairs.
#'
#' @param y matrix that stands for a sphere, S^d
#' @param i corresponds to the i-th observation index
#' @param rho param between 0 and 1 (not included)
#' @return spherical marginal Cauchy density probability given the parameters
#' @examples
#' y <- rotasym::r_unif_sphere(100, 2)
#' sc_cond_given_i(y, 1, 0.5)
#' sc_cond_given_i(y, 4, 0.9999)
sc_cond_given_i <- function(y, i, rho) {
  if (i < 1 || i > nrow(y)) {
    stop("The index i not valid, must be >= 1 and <= nrow(y)")
  }
  if (!rlang::is_scalar_atomic(rho)) {
    stop("Parameter rho must be an scalar")
  }
  # Calculate the probabilities of chosen a pair of elements where the i-th element is fixed
  return(sapply(1:nrow(y), jcondi_sc, y = y, i = i, rho))
}
