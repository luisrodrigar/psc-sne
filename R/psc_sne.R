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
  cosine_polysph <- cosine_polysph(x)
  rho_opt <- parLapply(cl, 1:n, function(i){
    optim(par = 0.5, 
          fn = function(rho) {
            res <- (to_perplexity(X = x, i = i, rho=rho, cosine_polysph = cosine_polysph) - perplexity)^2
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

simple_dspcauchy_hd <- function(x, i, j, rho, k, p) {
  ((1 + rho^2 - 2 * rho * t(x[i,,k]) %*% x[j,,k])^(-p))
}

P_ij_psc <- function(x, i, j, rho, p) {
  if(i == j)
    return(0)
  n <- nrow(x)
  r <- dim(x)[3]
  return(prod(sapply(1:r, FUN=simple_dspcauchy_hd, x=x, i=i, j=j, rho=rho, p=p)))
}

P_i_psc <- function(x, rho, p) {
  n <- nrow(x)
  r <- dim(x)[3]
  prob_is <- sapply(1:n, FUN=function(i){
    sapply(1:n, FUN=function(x, i, j){
      if(j!=i) {
        return(P_ij_psc(x, i, j, rho[i], p))
      } else {
        return(0)
      }
    }, x=x, i=i)
  }, simplify = 'array')
  return(rowSums(prob_is))
}


jcondi_psc <- function(x, i, j, rho, p, total_P=NULL) {
  if(is.null(total_P)) {
    total_P <- P_i_psc(x, rho, p)
  }
  return(P_ij_psc(x, i, j, rho[i], p)/total_P[i])
}


high_dimension_p <- function(X, rho_list, p) {
  X <- radial_projection_ps(X)
  n <- nrow(X)
  total_P <- P_i_psc(X, rho_list, p)
  jcondi <- function(i) {
    sapply(1:n, function(j) {
      jcondi_psc(X, i, j, rho_list, p, total_P)
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

kl_divergence_grad <- function(Y, i, rho, d, P) {
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


psc_sne <- function(X, p, d, rho_psc_list=NULL ,rho=0.5, perplexity=15, num_iteration=500, 
                    initial_momentum=0.5, final_momentum=0.8, eta=100, 
                    exageration=TRUE) {
  if(p+1!=ncol(X))
    stop("Error, the dimension of Y does not match with the value of p")
  if(d<1)
    stop("Error, d value must be greater or equal than 1")

  n <- nrow(X)
  if(is.null(rho_psc_list))
    rho_psc_list <- rho_optimize(X, perplexity)
  P_cond <- high_dimension_p(X, rho_list, p)
  P <- symmetric_probs(P_cond)
  
  total_iterations <- num_iteration + 2
  Y <- array(NA, c(n, d+1, total_iterations))
  Y[,,1] <- Y[,,2] <- r_unif_sphere(n, d+1)
  
  momentum = initial_momentum
  
  range_iterations <- seq_len(num_iteration) + 2
  for(i in range_iterations) {
    if(i>=0.75*num_iteration)
      momentum = final_momentum
    
    grad = t(sapply(1:n, kl_divergence_grad, Y=Y[,,i-1], rho=rho, d=d, P=P))
    ddir = - eta*grad
    
    Y_i <- Y[,,i-1] + ddir + momentum*(Y[,,i-1]-Y[,,i-2])
    Y[,,i] = radial_projection(Y_i)
    
    if(i %% 10 == 0) {
      Q <- low_dimension_Q(Y[,,i], d, rho)
      C = sum(P * log(P/Q), na.rm=TRUE)
      print(sprintf("Iteration %d: objective function value is %f", i, C))
    }
  }
  Y[,,total_iterations]
}

