library(DirStats)
library(Directional)
library(rotasym)
library(rlang)

## Polyspherical data generation

# r <- 3
# n <- 20
# d <- 2
gen_polysphere <- function(n, p, r) {
  polysphere <- array(NA, dim=c(n, p+1, r))
  for(k in seq_len(r)) {
    polysphere[,,k] <- r_unif_sphere(n, p+1)
  }
  polysphere
}

# polysphere <- gen_polysphere(n, p, r)

#######################################
###      Gaussian Distribution      ###
#######################################

# TODO: for the gaussian distribution which is implemented for the regular t-SNE

#####################################
### Spherical Cauchy Distribution ###
#####################################

## High-dimension 

simple_dspcauchy_hd <- function(x, i, j, rho, k, p) {
  if(i < 1 || i > nrow(x) || j < 1 || j > nrow(x))
    stop("The indexes (i and j) not valid, must be >= 1 and <= nrow(x)")
  if(k < 1 || k > dim(x)[3])
    stop("The 3rd dimensional index k not valid, must be >= 1 and <= r")
  if(!rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be an scalar")
  if(p+1 != ncol(x))
    stop("Parameter p is not valid, it must match with the dataset dimension")
  drop((1 + rho^2 - 2 * rho * t(x[i,,k]) %*% x[j,,k])^(-p))
}

P_ij_psc <- function(x, i, j, rho) {
  if(i < 1 || i > nrow(x) || j < 1 || j > nrow(x))
    stop("The indexes (i and j) not valid, must be >= 1 and <= nrow(x)")
  if(!rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be an scalar")
  if(i == j)
    return(0)
  return(prod(
    sapply(1:dim(x)[3], simple_dspcauchy_hd, x=x, i=i, j=j, rho=rho, p=ncol(x)-1)))
}

P_i_psc <- function(x, i, rho_list) {
  if(i < 1 || i > nrow(x))
    stop("The indexes i not valid, must be >= 1 and <= nrow(x)")
  if(rlang::is_scalar_atomic(rho_list))
    stop("Parameter rho must be a vector of dimension n")
  sapply(1:nrow(x), FUN=function(x, i, l, rho){
    if(l!=i) {
      return(P_ij_psc(x, i, l, rho))
    } else {
      return(0)
    }
  }, x=x, i=i, rho=rho_list[i])
}

P_total_psc <- function(x, rho_list) {
  if(rlang::is_scalar_atomic(rho_list))
    stop("Parameter rho must be a vector of dimension n")
  return(rowSums(sapply(1:nrow(x), 
                        FUN=P_i_psc, x=x, rho=rho_list, 
                        simplify = 'array')))
}

jcondi_psc <- function(x, i, j, rho_list, total_P_i) {
  if(i < 1 || i > nrow(x) || j < 1 || j > nrow(x))
    stop("The indexes (i and j) not valid, must be >= 1 and <= nrow(x)")
  if(rlang::is_scalar_atomic(rho_list))
    stop("Parameter rho must be a vector of dimension n")
  if(!rlang::is_scalar_atomic(total_P_i))
    stop("Parameter total_P_i must be an scalar")
  if(i==j)
    return(0)
  return(P_ij_psc(x, i, j, rho_list[i])/total_P_i)
}

psc_cond_given_i <- function(x, i, rho_list, total_P_i=NULL) {
  if(i < 1 || i > nrow(x))
    stop("The index i not valid, must be >= 1 and <= nrow(x)")
  if(rlang::is_scalar_atomic(rho_list))
    stop("Parameter rho must be a vector of dimension n")
  if(!rlang::is_scalar_atomic(total_P_i))
    stop("Parameter total_P_i must be an scalar")
  if(is.null(total_P_i)) {
    total_P_i <- sum(P_i_psc(x, i, rho_list))
  }
  return(sapply(1:nrow(x), jcondi_psc, x=x, i=i, rho_list=rho_list, 
                total_P_i=total_P_i))
}

## Low-dimension

simple_dspcauchy_ld <- function(y, i, j, rho, d) {
  if(i < 1 || i > nrow(y) || j < 1 || j > nrow(y))
    stop("The indexes (i and j) not valid, must be >= 1 and <= nrow(y)")
  if(!rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be an scalar")
  if(d+1 != ncol(y))
    stop("Parameter p is not valid, it must match with the dataset dimension")
  drop((1 + rho^2 - 2 * rho * t(y[i,]) %*% y[j,])^(-d))
}

Q_i_sc <- function(y, i, rho) {
  if(i < 1 || i > nrow(y))
    stop("The index i not valid, must be >= 1 and <= nrow(y)")
  if(!rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be an scalar")
  sum(sapply((1:nrow(y))[-i], simple_dspcauchy_ld, y=y, i=i, rho=rho, d=ncol(y)-1))
}

jcondi_sc <- function(y, i, j, rho, total_Q_i=NULL) {
  if(i < 1 || i > nrow(y))
    stop("The index i not valid, must be >= 1 and <= nrow(y)")
  if(!rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be an scalar")
  if(is.null(total_Q_i))
    total_Q_i = Q_i_sc(y, i, rho)
  if(i==j)
    return(0)
  return(simple_dspcauchy_ld(y, i, j, rho, ncol(y)-1)/total_Q_i)
}

sc_cond_given_i <- function(y, i, rho) {
  if(i < 1 || i > nrow(y))
    stop("The index i not valid, must be >= 1 and <= nrow(y)")
  if(!rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be an scalar")
  return(sapply(1:nrow(y), jcondi_sc, y=y, i=i, rho))
}




