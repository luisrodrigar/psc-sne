library(DirStats)
library(Directional)
library(rotasym)

## Polyspherical data generation

# r <- 3
# n <- 20
# d <- 2
gen_polysphere <- function(n, d, r) {
  p <- (d+1)
  polysphere <- array(NA, dim=c(n, p, r))
  for(k in seq_len(r)) {
    polysphere[,,k] <- r_unif_sphere(n, p)
  }
  polysphere
}

# polysphere <- gen_polysphere(n, d, r)

#######################################
###      Gaussian Distribution      ###
#######################################

# TODO: for the gaussian distribution which is implemented for the regular t-SNE

#####################################
### Spherical Cauchy Distribution ###
#####################################

## High-dimension 

simple_dspcauchy <- function(x, i, j, rho, k, p) {
  ((1 + rho^2 - 2 * rho * t(x[i,,k]) %*% x[j,,k])^(-p))
}

P_ij_psc <- function(x, i, j, rho, p) {
  if(i == j)
    return(0)
  n <- nrow(x)
  r <- dim(x)[3]
  return(prod(sapply(1:r, FUN=simple_dspcauchy, x=x, i=i, j=j, rho=rho, p=p)))
}

simple_dspcauchy_pov <- function(x, i, j, rho, p) {
  ((1 + rho^2 - 2 * rho * t(x[i,]) %*% x[j,])^(-p))
}

P_ij_psc_pov <- function(x, i, j, rho, p, r) {
  if(i == j)
    return(0)
  n <- nrow(x)
  r <- dim(x)[3]
  return(prod(sapply(1:r, FUN=simple_dspcauchy, x=x, i=i, j=j, rho=rho, p=p)))
}

P_i_psc <- function(x, rho, d) {
  n <- nrow(x)
  r <- dim(x)[3]
  prob_is <- sapply(1:n, FUN=function(i){
    sapply(1:n, FUN=function(x, i, j){
      if(j!=i) {
        return(p_ij_sc(x, i, j, rho[i], d))
      } else {
        return(0)
      }
    }, x=x, i=i)
  }, simplify = 'array')
  return(rowSums(prob_is))
}

jcondi_sc <- function(x, i, j, rho, d, total_p=NULL) {
  if(is.null(total_p)) {
    total_p <- p_i_sc(x, rho, d)
  }
  return(p_ij_sc(x, i, j, rho[i], d)/total_p[i])
}


## Low-dimension

simple_dspcauchy_ld <- function(x, i, j, rho, d) {
  ((1 + rho^2 - 2 * rho * t(x[i,]) %*% x[j,])^(-d))
}

jcondi_spcauchy_ld <- function(Y, i, rho, d) {
  prob_i <- function(Y, i, l, rho, d) {
    simple_dspcauchy_ld(x=Y, i=i, j=l, rho, d)
  }
  prob_ij <- sapply(1:n, simple_dspcauchy_ld, x=Y, i=i, rho, d)
  (prob_ij / sum(sapply(1:n, prob_i, Y=Y, i=i, rho=rho, d=d)))
}




