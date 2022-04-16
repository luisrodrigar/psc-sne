library(DirStats)
library(Directional)
library(rotasym)
library(parallel)
library(doParallel)

## Polyspherical data generation

k <- 2
n <- 100
d <- 2
gen_polysphere <- function(n, d, r) {
  p <- (d+1)
  polysphere <- array(NA, dim=c(n, p, r))
  for(k in seq_len(r)) {
    polysphere[,,k] <- r_unif_sphere(n, p)
  }
  polysphere
}

#######################################
###      Gaussian Distribution      ###
#######################################

# TODO: for the gaussian distribution which is implemented for the regular t-SNE

#####################################
### Spherical Cauchy Distribution ###
#####################################

## High-dimension 

simple_dspcauchy <- function(x, l, i, rho, k, d) {
  ((1 + rho^2 - 2 * rho * t(x[l,,k]) %*% x[i,,k])^(-d))
}

prob_i_spcauchy <- function(x, i, rho, d){
  r <<- dim(x)[3]
  n <<- nrow(x)
  probi_k_spcauchy <- function(j) {
    prod(sapply(X=seq_len(r), FUN=simple_dspcauchy, x=x, i=i, l=j, rho=rho, d=d))
  }
  sum(sapply(seq_len(n)[-i], probi_k_spcauchy))
}

jcondi_spcauchy <- function(x, i, j, rho, d, prob_is=NULL) {
  r <- dim(x)[3]
  if(!is.null(prob_is) && !is.null(prob_is[i])) {
    prob_i <<- prob_is[i]
  } else {
    prob_i <<- prob_i_spcauchy(x, i, rho, d)
  }
  (prod(sapply(1:r, simple_dspcauchy, x=x, i=i, l=j, rho=rho, d=d)) / prob_i)
}

## Low-dimension

simple_dspcauchy_ld <- function(x, l, i, rho, d) {
  ((1 + rho^2 - 2 * rho * t(x[l,]) %*% x[i,])^(-d))
}

jcondi_spcauchy_ld <- function(Y, i, rho, d) {
  prob_i <- function(Y, i, l, rho, d) {
    simple_dspcauchy_ld(x=Y, l=l, i=i, rho, d)
  }
  prob_ij <- sapply(1:n, simple_dspcauchy_ld, x=Y, i=i, rho, d)
  (prob_ij / sum(sapply(1:n, prob_i, Y=Y, i=i, rho=rho, d=d)))
}





