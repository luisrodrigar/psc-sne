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

diag_3d <- function(x, k, val) {
  diag(x[,,k]) <- val
  return(x[,,k])
}

## cosine similarity

cosine_polysph <- function(X) {
  library(lsa)
  
  r <- dim(X)[3]
  cosine_sphere_ith <- function(k){
    cos_sim <- cosine(t(X[,,k]))
  }
  sapply(1:r, cosine_sphere_ith, simplify = 'array')
}

## optimal version

high_dimension <- function(x, rho, cosine_polysphere=NULL) {
  n <- nrow(x)
  d <- ncol(x)-1
  r <- dim(x)[3]
  cos_sim_pol <- cosine_polysph(x)
  
  P <- sweep(cos_sim_pol, MARGIN=1, STATS=(-2*rho), FUN="*", 
             check.margin=FALSE)
  P <- sweep(P, MARGIN=1, STATS=(rho^2), FUN="+", 
             check.margin=FALSE)
  P <- 1/(1+P)^d
  
  Paux <- P
  P <- sapply(1:r, FUN=diag_3d, x=Paux, val=0, simplify = 'array')
  P_i_r <- apply(P, MARGIN=c(1,2), prod)
  Pi <- rowSums(P_i_r)
  P_ij = sweep(x=P_i_r, MARGIN=c(1,2), STATS=Pi, FUN="/", 
               check.margin=FALSE)
  return(P_ij)
}

## inefficient version

high_dimension_p <- function(X, rho_list, d) {
  n <- nrow(X)
  total_p <- p_i_sc(X, rho_list, d)
  jcondi <- function(i) {
    sapply(1:n, function(j) {
      jcondi_sc(X, i, j, rho_list, d, total_p)
    })
  }
  P <- sapply(1:n, jcondi)
  return(t(P))
}
