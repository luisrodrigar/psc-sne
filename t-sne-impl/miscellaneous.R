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

jcondi_spcauchy <- function(x, i, j, rho, d) {
  r <- dim(x)[3]
  (prod(sapply(1:r, simple_dspcauchy, x=x, i=i, l=j, rho=rho, d=d)) / 
      prob_i_spcauchy(x, i, rho, d))
}

to_perplexity <- function(X, i, rho, d=2) {
  entropy <- function(j) {
    jcondi_value <- jcondi_spcauchy(X, i, j, rho, d)
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

rho_opt <- parLapply(cl, 1:n, function(i){
  print(paste("The ", i, "-th observation", sep = ""))
  optim(par = 0.5, 
        fn = function(rho) {
          res <- (to_perplexity(X = polysphere, i = i, rho=rho) - perplex)^2
          ifelse(is.finite(res), res, 1e6)
        },
        method="L-BFGS-B", lower=0, upper=.9999)$par
})
rho_opt <- simplify2array(rho_opt)
stopCluster(cl) # 3min-ish in parallel
