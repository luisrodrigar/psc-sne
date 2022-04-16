###############################################
##        Poly-spherical Cauchy HD           ##
## High-dimension neighborhood probabilities ##
###############################################

# utils

symmetric_probs <- function(P) {
  n <- nrow(P)
  P = (P + t(P)) / (2*n)
  return(P)
}

diag_3d <- function(array3d, k, value) {
  diag(array3d[,,k]) <- value
  array3d[,,k]
}

## cosine similarity

library(lsa)

cosine_polysph <- function(polysphere) {
  r <- dim(polysphere)[3]
  cosine_sphere_ith <- function(k){
    cos_sim <- cosine(t(polysphere[,,k]))
  }
  sapply(1:r, cosine_sphere_ith, simplify = 'array')
}

## optimal version

high_dimension <- function(x, rho) {
  n <- nrow(x)
  d <- ncol(x)-1
  r <- dim(x)[3]
  cos_sim_pol <- cosine_polysph(x)
  
  P <- sweep(cos_sim, MARGIN=1, STATS=(-2*rho), FUN="*")
  P <- sweep(P, MARGIN=1, STATS=(rho^2), FUN="+")
  P <- 1/(1+P)^d
  
  P <- sapply(1:r, FUN=diag_3d, array3d=P, value=0, simplify = 'array')
  P_i_r <- apply(P, MARGIN=c(1,2), prod)
  Pi <- rowSums(P_i_r)
  P_ij = sweep(x=P_i_r, MARGIN=c(1,2), STATS=Pi, FUN="/")
  return(P_ij)
}

## inefficient version

high_dimension_P <- function(X, d, rho) {
  n <- nrow(X)
  jcondi <- function(i) {
    prob_is <- sapply(seq_len(n), prob_i_spcauchy, x=X, rho=rho[i], d=d)
    sapply(1:n, function(j) {
      jcondi_spcauchy(X, i, j, rho[i], d, prob_is)
    })
  }
  P <- sapply(1:n, jcondi)
  diag(P) <- 0
  return(t(P))
}
P <- high_dimension_P(polysphere, 2, rho_opt)

Pij <- symmetric_probs(P)