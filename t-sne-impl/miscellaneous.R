###############################
###############################
##                           ##
## Perplexity Generalization ##
##                           ##
###############################
###############################

#######################################
###      Gaussian Distribution      ###
#######################################

# TODO: for the gaussian distribution which is implemented for the regular t-SNE

#####################################
### Spherical Cauchy Distribution ###
#####################################

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

to_perplexity <- function(X, i, rho, d=2, prob_is=NULL) {
  entropy <- function(j) {
    jcondi_value <- jcondi_spcauchy(X, i, j, rho, d, prob_is)
    (jcondi_value * log2(jcondi_value))
  }
  return(2^(-1*sum(sapply(seq_len(nrow(X))[-i], entropy))))
}

# Polyspherical data generation

library(DirStats)
library(Directional)

k <- 2
n <- 100
polysphere <- array(NA, dim=c(n,3, k))
for(i in seq_len(k)) {
  th <- sample(seq(0, pi/2, l=200), size=n)
  ph <- sample(seq(0, 2*pi, l=200), size=n)
  polysphere[,,i] <- to_sph(th, ph)
}

library(parallel)
library(doParallel)

perplex <- 25
num_cores <- detectCores()-1
cl <- makeForkCluster(num_cores)
start_time <- Sys.time()
rho_opt <- parLapply(cl, 1:n, function(i){
  optim(par = 0.5, 
        fn = function(rho) {
          prob_is <- sapply(seq_len(n), prob_i_spcauchy, x=polysphere, rho=rho, d=2)
          res <- (to_perplexity(X = polysphere, i = i, rho=rho, prob_is=prob_is) - perplex)^2
          ifelse(is.finite(res), res, 1e6)
        },
        method="L-BFGS-B", lower=0, upper=.9999)$par
})
end_time <- Sys.time()
print(end_time-start_time)
rho_opt <- simplify2array(rho_opt)
stopCluster(cl) # Time difference of 3.92566 mins

## P

high_dimension_P <- function(X, d, rho) {
  n <- nrow(X)
  jcondi <- function(i) {
    prob_is <- sapply(seq_len(n), prob_i_spcauchy, x=polysphere, rho=rho[i], d=d)
    sapply(1:n, function(j) {
      jcondi_spcauchy(X, i, j, rho[i], d, prob_is)
    })
  }
  P <- sapply(1:n, jcondi)
  diag(P) <- 0
  return(t(P))
}
P <- high_dimension_P(polysphere, 2, rho_opt)
