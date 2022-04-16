##############################################
##        Poly-spherical Cauchy HD          ##
## Low-dimension neighborhood probabilities ##
##############################################

p <- d+1
total_iterations <- 250
Y <- array(NA, c(n, p, total_iterations))
Y[,,1] <- Y[,,2] <- r_unif_sphere(n, p)

low_dimension_Q <- function(Y, d, rho, total_iterations) {
  cos_simil <- cosine(t(Y))
  Q <- 1/(1+rho^2-2*rho*cos_simil)^d
  diag(Q) <- 0
  Qi <- rowSums(Q)
  Q_ij = sweep(x=Q, MARGIN=1, STATS=Qi, FUN="/")
  return(Q_ij)
}

Q <- low_dimension_Q(polysphere, 2, 0.5, 250)
Qij <- symmetric_probs(Q)