library(lsa)
library(parallel)
library(doParallel)

cos_sim_ij <- function(i, j) {
  (t(i) %*% j)[1]
}

## High-dimensional space

diag_3d <- function(x, k, val) {
  diag(x[,,k]) <- val
  return(x[,,k])
}

### cosine similarity

cosine_polysph <- function(X) {
  r <- dim(X)[3]
  cosine_sphere_ith <- function(k){
    cos_sim <- cosine(t(X[,,k]))
  }
  sapply(1:r, cosine_sphere_ith, simplify = 'array')
}

### Radial projection

radial_projection <- function(X) {
  t(sapply(1:nrow(X), function(i){
    X[i,]/norm(X[i,], type="2")
  }))
}

radial_projection_ps <- function(X) {
  r = dim(X)[3]
  sapply(1:r, function(k) radial_projection(X[,,k]), simplify = 'array')
}

### Perplexity

rho_optimize <- function(x, perplexity) {
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

simple_dspcauchy_hd <- function(x, i, j, rho, k, d) {
  ((1 + rho^2 - 2 * rho * t(x[i,,k]) %*% x[j,,k])^(-d))
}

p_ij_sc <- function(x, i, j, rho, d) {
  if(i == j)
    return(0)
  n <- nrow(x)
  r <- dim(x)[3]
  return(prod(sapply(1:r, FUN=simple_dspcauchy_hd, x=x, i=i, j=j, rho=rho, d=d)))
}

p_i_sc <- function(x, rho, d) {
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


high_dimension_p <- function(X, rho_list, d) {
  X <- radial_projection_ps(X)
  n <- nrow(X)
  total_p <- p_i_sc(X, rho_list, d)
  jcondi <- function(i) {
    sapply(1:n, function(j) {
      jcondi_sc(X, i, j, rho_list, d, total_p)
    })
  }
  P <- sapply(1:n, jcondi)
  return(t(P))
}

symmetric_probs <- function(P) {
  n <- nrow(P)
  P = (P + t(P)) / (2*n)
  return(P)
}

low_dimension_Q <- function(Y, d, rho) {
  Z <- radial_projection(Y)
  cos_simil <- cosine(t(Z))
  Q <- (1+rho^2-2*rho*cos_simil)^(-d)
  diag(Q) <- 0
  Qi <- sum(Q)
  Q_ij = Q/Qi
  return(Q_ij)
}

kl_divergence_grad <- function(Y, i,  rho, d, P) {
  if(i < 1 || i > nrow(Y))
    stop(paste("Error, i value not allowed. Positive values greater tha 0 and",
         "smaller or equal than the total number of observations."))
  if(d<1)
    stop("Error, d value not allowed. Positive values greater or equal than 1.")
  if(rho < 0 || rho >= 1)
    stop(paste("Error, rho value not allowed, values in between 0 and 1,", 
               "the last one not included."))
  if(nrow(Y)!=nrow(P))
    stop("Error, the num of observations does not match with the matrix P")
  if(d+1!=ncol(Y))
    stop("Error, the columns of Y does not match with the value of d")
  
  Z <- radial_projection(Y)
  Q <- low_dimension_Q(Z, d, rho)
  n_minus_i <- (1:nrow(Y))[-i]
  (4*d*rho*colSums(t(
    sapply(n_minus_i, function(j) {
      s_ij <- cos_sim_ij(Z[i,], Z[j,])
      Q_minus_P <- (Q[i,j]-P[i,j])
      div <- (1+rho^2-2*rho*s_ij)
      return(Z[j,]/div*Q_minus_P)
    }))))
}