
## proposal

entropy_func <- function(X, i, var) {
  beta = 1/(2*var)
  P_i <- rowSums(exp(-(t(t(X) - X[i, ])^2) * beta))
  return (log(sum(P_i)) + (beta * sum((t(t(X) - X[i, ])^2) * P_i) / sum(P_i)))
}

calc_perplexity <- function(X, i, var) {
  exp(entropy_func(X, i, var))
}

diff_perpl_optim <- function(X, i, perplex, var) {
  res <- (calc_perplexity(X, i, var) - perplex)^2
  if(is.finite(res))
    return(res)
  else
    return(1e6)
}

perpl_to_var <- function(X, perplexity) {
  sapply(1:n, function(i) {
    optim_func <- partial(diff_perpl_optim, X=X, i=i, perplex = perplexity)
    optim(par = 0.75, fn = optim_func, method = "L-BFGS-B", lower = 0.1)$par
  })
}

## tfg

calc_perplexity <-function(X, i, sigma2) {
  
  # Dimensions
  n <- nrow(X)
  p <- ncol(X)
  
  P_i_cond_j= softmax(- rowSums( t(t(X) - X[i, ])^2) / (2 * sigma2))
  H_i <- -sum(P_i_cond_j * log2(P_i_cond_j))
  
  #returns the perplexity
  return(2^H_i)
  
}

perplex = 30
X<- data

sapply(1:n, function(i) {
  
  optim(par = 0.75, 
        fn = function(s2) {
          res <- (calc_perplexity(X = X, i = i, sigma2 = s2) - perplex)^2
          ifelse(is.finite(res), res, 1e6)
        },
        method = "L-BFGS-B", lower = 0.1)$par
  
})