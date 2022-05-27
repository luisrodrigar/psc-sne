library(parallel)
library(doParallel)
library(lsa)

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
  sapply(1:r, function(x, k) cosine(t(x[,,k])), x=X, simplify = 'array')
}

### Radial projection

radial_projection <- function(X) {
  t(sapply(1:nrow(X), function(i){
    X[i,]/norm(X[i,], type="2")
  }))
}

radial_projection_ps <- function(X) {
  r = dim(X)[3]
  sapply(1:r, function(x, k) radial_projection(x[,,k]), x=X, simplify = 'array')
}

### Perplexity

rho_optimize <- function(x, perplexity, cosine_polysph=NULL, num_cores = 2) {
  n <- nrow(x)
  cl <- makeForkCluster(num_cores)
  start_time <- Sys.time()
  if(is.null(cosine_polysph))
    cosine_polysph <- cosine_polysph(x)
  rho_opt <- parLapply(cl, 1:n, function(i){
    stats::optim(par = 0.5, 
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
  if(is.null(total_P))
    total_P <- P_i_psc(x, rho, p)
  return(P_ij_psc(x, i, j, rho[i], p)/total_P[i])
}


high_dimension <- function(x, rho, cosine_polysphere=NULL) {
  n <- nrow(x)
  d <- ncol(x)-1
  r <- dim(x)[3]
  if(is.null(cosine_polysphere))
    cosine_polysphere <- cosine_polysph(x)
  
  P <- sweep(cosine_polysphere, MARGIN=1, STATS=(-2*rho), FUN="*", 
             check.margin=FALSE)
  P <- sweep(P, MARGIN=1, STATS=(rho^2), FUN="+", 
             check.margin=FALSE)
  P <- 1/(1+P)^d
  
  Paux <- P
  P <- sapply(1:r, FUN=diag_3d, x=Paux, val=0, simplify = 'array')
  P_i_r <- apply(P, MARGIN=c(1,2), prod)
  Pi <- rowSums(P_i_r)
  P_ij = sweep(x=P_i_r, MARGIN=c(1,2), STATS=Pi, FUN="/", 
               check.margin=FALSE)
  return(P_ij)
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

kl_divergence_grad <- function(Y, i, rho, d, P, cos_sim=NULL) {
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
  if(is.null(cos_sim))
    cos_sim = cosine(t(Z))
  Q <- low_dimension_Q(Z, d, rho)
  n_minus_i <- (1:nrow(Y))[-i]
  (4*d*rho*colSums(t(
    sapply(n_minus_i, function(j) {
      return(Z[j,]/(1+rho^2-2*rho*cos_sim[i,j])*(Q[i,j]-P[i,j]))
    }))))
}


psc_sne <- function(X, d, rho_psc_list=NULL ,rho=0.5, perplexity=15, num_iteration=200, 
                    initial_momentum=0.5, final_momentum=0.8, eta=100, 
                    exageration=TRUE, colors=NULL) {
  if(d<1)
    stop("Error, d value must be greater or equal than 1")

  n <- nrow(X)
  p <- ncol(X)-1
  
  cosine_sim_polysphere <- cosine_polysph(X)
  if(is.null(rho_psc_list))
    rho_psc_list <- rho_optimize(X, perplexity, 
                                 cosine_polysph = cosine_sim_polysphere)
  P_cond <- high_dimension(x=X, rho=rho_psc_list, cosine_polysphere = cosine_sim_polysphere)
  P <- symmetric_probs(P_cond)
  
  total_iterations <- num_iteration + 2
  Y <- array(NA, c(n, d+1, total_iterations))
  Y[,,1] <- Y[,,2] <- r_unif_sphere(n, d+1)
  
  momentum = initial_momentum
  max_tol = 1e-3
  
  range_iterations <- seq_len(num_iteration) + 2
  for(i in range_iterations) {
    print(sprintf("Iteration %d", i))
    if(i>=0.75*num_iteration)
      momentum = final_momentum
    
    cos_sim <- cosine(t(Y[,,i-1]))
    grad = t(simplify2array(mclapply(mc.cores = detectCores()-1, 1:n, 
                      kl_divergence_grad, Y=Y[,,i-1], rho=rho, d=d, P=P, 
                      cos_sim = cos_sim)))
    ddir = - eta*grad
    
    Y_i <- Y[,,i-1] + ddir + momentum*(Y[,,i-1]-Y[,,i-2])
    Y[,,i] = radial_projection(Y_i)
    
    if(i == 3 || (i-2) %% 5 == 0) {
      Q <- low_dimension_Q(Y[,,i], d, rho)
      C = sum(P * log(P/Q), na.rm=TRUE)
      print(sprintf("Iteration %d: objective function value is %f", i, C))
      if(d==1) {
        Y_rad <- DirStats::to_rad(Y[,,i])
        r <- 1
        theta <<- Y_rad
        if(is.null(colors))
          colors <- rep(1, n)
        plot(r*sin(theta),
             r*cos(theta),
             col=colors,
             xlim=c(-max(r),max(r)),
             ylim=c(-max(r),max(r)))
        
        polygon(max(r)*sin(seq(0,2*pi,length.out=100)),max(r)*cos(seq(0,2*pi,length.out=100)))
      }
    }
    tol <- norm(Y[,,i]-Y[,,i-1], "2")
    if(tol < max_tol)
      break
  }
  Y[,,total_iterations]
}

