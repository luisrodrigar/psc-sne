###############################################
##        Poly-spherical Cauchy HD           ##
## High-dimension neighborhood probabilities ##
###############################################

library(lsa)

# utils

symmetric_probs <- function(P) {
  n <- nrow(P)
  P <- (P + t(P)) / (2 * n)
  return(P)
}

diag_3d <- function(x, k, val) {
  if (k < 1 || k > dim(x)[3]) {
    stop("The 3rd dimensional index k not valid, must be >= 1 and <= r")
  }
  diag(x[, , k]) <- val
  return(x[, , k])
}

## cosine similarity

cosine_polysph <- function(x) {
  r <- dim(x)[3]
  cosine_sphere_ith <- function(k) {
    cos_sim <- cosine(t(x[, , k]))
  }
  sapply(1:r, cosine_sphere_ith, simplify = "array")
}

## optimal version

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
  n <- nrow(x)
  d <- ncol(x) - 1
  r <- dim(x)[3]
  if (is.null(cos_sim_pol)) {
    cos_sim_pol <- cosine_polysph(x)
  }

  P <- sweep(cos_sim_pol,
    MARGIN = 1, STATS = (-2 * rho_list), FUN = "*",
    check.margin = FALSE
  )
  P <- sweep(P,
    MARGIN = 1, STATS = (rho_list^2), FUN = "+",
    check.margin = FALSE
  )
  P <- 1 / (1 + P)^d

  Paux <- P
  P <- sapply(1:r, FUN = diag_3d, x = Paux, val = 0, simplify = "array")
  P_i_r <- apply(P, MARGIN = c(1, 2), prod)
  Pi <- rowSums(P_i_r)
  P_ij <- sweep(
    x = P_i_r, MARGIN = c(1, 2), STATS = Pi, FUN = "/",
    check.margin = FALSE
  )
  return(P_ij)
}

## inefficient version

high_dimension_p <- function(x, rho_list) {
  if (!rlang::is_vector(rho_list)) {
    stop("Parameter rho_list must be a vector")
  }
  if (length(rho_list) != nrow(x)) {
    stop("Parameter rho_list not valid, size must be equal to nrow(x)")
  }
  total_p <- P_total_psc(x, rho_list)
  jcondi <- function(i) {
    sapply(1:nrow(x), function(j) {
      jcondi_psc(x, i, j, rho_list, total_p[i])
    })
  }
  return(t(sapply(1:nrow(x), jcondi)))
}
