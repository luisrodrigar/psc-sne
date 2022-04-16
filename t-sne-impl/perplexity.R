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

perplex <- 30
num_cores <- detectCores()-1
cl <- makeForkCluster(num_cores)
start_time <- Sys.time()
rho_opt <- parLapply(cl, 1:n, function(i){
  optim(par = 0.5, 
        fn = function(rho) {
          prob_is <- sapply(seq_len(n), prob_i_spcauchy, x=polysphere, rho=rho)
          res <- (to_perplexity_P(X = polysphere, i = i, rho=rho, prob_is=prob_is) - perplex)^2
          ifelse(is.finite(res), res, 1e6)
        },
        method="L-BFGS-B", lower=0, upper=.9999)$par
})
end_time <- Sys.time()
print(end_time-start_time)
rho_opt <- simplify2array(rho_opt)
stopCluster(cl) # Time difference of 21.67933 mins

perplex <- 30
num_cores <- detectCores()-1
cl <- makeForkCluster(num_cores)
start_time <- Sys.time()
rho_opt <- parLapply(cl, 1:n, function(i){
  optim(par = 0.5, 
        fn = function(rho) {
          res <- (to_perplexity(X = polysphere, i = i, rho=rho) - perplex)^2
          ifelse(is.finite(res), res, 1e6)
        },
        method="L-BFGS-B", lower=0, upper=.9999)$par
})
end_time <- Sys.time()
print(end_time-start_time)
rho_opt <- simplify2array(rho_opt)
stopCluster(cl) # Time difference of 5.836471 mins


to_perp <- function(X, rho) {
  P <- high_dimension(X, rho)
  op <- P*log2(P)
  diag(op) <- 0
  return(2^(-rowSums(op)))
}

n <- nrow(polysphere)
perps <- rep(30, n)
start_time <- Sys.time()
optim(par = rep(0.5, n), 
      fn = function(rho) {
        print(rho)
        res <- norm((to_perp(polysphere, rho) - perps)^2, type="2")
        print(res)
        ifelse(is.finite(res), res, 1e6)
      },
      method="L-BFGS-B", lower=rep(0, n), upper=rep(.9999, n))$par
end_time <- Sys.time()
print(end_time-start_time) # Time difference of 4.311193 mins

