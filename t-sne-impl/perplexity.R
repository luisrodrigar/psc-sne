###########################################
##      spherical Cauchy Perplexity      ##
###########################################

to_perplexity_P <- function(X, i, rho, prob_is=NULL) {
  n <- nrow(X)
  d <- (ncol(X)-1)
  entropy <- function(j) {
    jcondi_value <- jcondi_spcauchy(X, i, j, rho, prob_is)
    (jcondi_value * log2(jcondi_value))
  }
  return(2^(-1*sum(sapply(seq_len(n)[-i], entropy))))
}

to_perplexity <- function(X, i, rho) {
  n <- nrow(X)
  d <- (ncol(X)-1)
  Picondj <- high_dimension(X, rep(rho, n))
  entropy <- function(j) {
    (Picondj[i,j] * log2(Picondj[i, j]))
  }
  return(2^(-1*sum(sapply(seq_len(n)[-i], entropy))))
}

to_perp <- function(X, rho) {
  P <- high_dimension(X, rho)
  op <- P*log2(P)
  diag(op) <- 0
  return(2^(-rowSums(op)))
}

# inefficient process to calculate the rho

## scalar calculus

rho_optim_inefficient <- function(x, perplexity) {
  n <- nrow(x)
  num_cores <- detectCores()-1
  cl <- makeForkCluster(num_cores)
  start_time <- Sys.time()
  rho_opt <- parLapply(cl, 1:n, function(i){
    optim(par = 0.5, 
          fn = function(rho) {
            prob_is <- sapply(seq_len(n), prob_i_spcauchy, x=polysphere, rho=rho)
            res <- (to_perplexity_P(X = polysphere, i = i, rho=rho, prob_is=prob_is) - perplexity)^2
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

rho_optim_inefficient <- function(x, perplexity) {
  n <- nrow(x)
  perps <- rep(perplexity, n)
  start_time <- Sys.time()
  rho_opt <- optim(par = rep(0.5, n), 
                   fn = function(rho) {
                     res <- norm(x=((to_perp(polysphere, rho) - perps)^2), type="2")
                     ifelse(is.finite(res), res, 1e6)
                   },
                   method="L-BFGS-B", lower=rep(0, n), upper=rep(.9999, n))$par
  end_time <- Sys.time()
  print(end_time-start_time) 
  return(rho_opt)
}

## efficient ways to calculate the perplexity

## parallel and matrix calculus 

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



