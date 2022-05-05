##############################################
##        Poly-spherical Cauchy HD          ##
## Low-dimension neighborhood probabilities ##
##############################################

library(lsa)

low_dimension_Q <- function(Y, d, rho) {
  cos_simil <- cosine(t(Y))
  Q <- (1+rho^2-2*rho*cos_simil)^(-d)
  diag(Q) <- 0
  Qi <- sum(Q)
  Q_ij = Q/Qi
  return(Q_ij)
}