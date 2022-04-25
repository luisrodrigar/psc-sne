###########################################
##      spherical Cauchy Perplexity      ##
###########################################
library(parallel)

to_perplexity_P <- function(X, i, rho) {
  n <- nrow(X)
  d <- (ncol(X)-1)
  rho_list <- rep(rho, n)
  total_p <- p_i_sc(X, rho_list, d)
  Picondj <- sapply(seq_len(n), jcondi_sc, x=X, i=i, rho=rho_list, d=d, 
                    total_p=total_p)
  entropy <- function(j) {
    return(Picondj[j] * log2(Picondj[j]))
  }
  return(2^(-sum(sapply(seq_len(n)[-i], entropy))))
}

to_perplexity <- function(X, i, rho) {
  n <- nrow(X)
  d <- (ncol(X)-1)
  cosine_polysph <- cosine_polysph(X)
  Picondj <- high_dimension(X, rep(rho, n), cosine_polysph)
  entropy <- function(j) {
    (Picondj[i,j] * log2(Picondj[i,j]))
  }
  return(2^(-sum(sapply(seq_len(n)[-i], entropy))))
}

to_perp <- function(X, rho) {
  cosine_polysph <- cosine_polysph(X)
  P <- high_dimension(X, rho, cosine_polysph)
  op <- P*log2(P)
  diag(op) <- 0
  return(2^(-rowSums(op)))
}

# inefficient process to calculate the rho

## scalar calculus
## Time difference of 25.88906 mins
## 50 rows 25 spheres
rho_optim_inefficient <- function(x, perplexity) {
  n <- nrow(x)
  d <- (ncol(x)-1)
  num_cores <- detectCores()-1
  cl <- makeForkCluster(num_cores, outfile="log.txt")
  start_time <- Sys.time()
  rho_opt <- parLapply(cl, 1:n, function(i){
    print(i)
    optim(par = 0.5, 
          fn = function(rho) {
            total_p <- p_i_sc(x, rep(rho, n), d)
            res <- (to_perplexity_P(x, i, rho) - perplexity)^2
            print(res)
            ifelse(is.finite(res), res, 1e6)
          },
          method="L-BFGS-B", lower=0, upper=.9999)$par
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  rho_opt <- simplify2array(rho_opt)
  stopCluster(cl)
  return(rho_opt)
}

## matricidal calculus, cost function
## Time difference of 14.22428 mins
## 50 rows 25 spheres

rho_optim_ineff <- function(x, perplexity) {
  n <- nrow(x)
  num_cores <- detectCores()-1
  cl <- makeForkCluster(num_cores, outfile="log.txt")
  start_time <- Sys.time()
  rho_opt <- parLapply(cl, 1:n, function(i){
    print(paste("Observation :", i))
    optim(par = 0.5, 
          fn = function(rho) {
            res <- (to_perp(polysphere, rep(rho,n))[i] - perplexity)^2
            print(res)
            ifelse(is.finite(res), res, 1e6)
            },
          method="L-BFGS-B", lower=0, upper=.9999)$par
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  rho_opt <- simplify2array(rho_opt)
  stopCluster(cl)
  return(rho_opt)
}

## efficient ways to calculate the perplexity

## parallel and matrix calculus 
## Time difference of 14.56007 mins
## 50 rows 25 spheres

rho_optimize <- function(x, perplexity) {
  library(parallel)
  library(doParallel)
  library(lsa)
  n <- nrow(x)
  num_cores <- detectCores()-1
  cl <- makeForkCluster(num_cores)
  start_time <- Sys.time()
  rho_opt <- parLapply(cl, 1:n, function(i){
    optim(par = 0.5, 
          fn = function(rho) {
            res <- (to_perplexity(X = x, i = i, rho=rho) - perplexity)^2
            ifelse(is.finite(res), res, 1e6)
          },
          method="L-BFGS-B", lower=0, upper=.9999)$par
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  rho_opt <- simplify2array(rho_opt)
  stopCluster(cl) # Time difference of 5.836471 mins
  return(rho_opt)
}



