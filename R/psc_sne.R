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

low_dimension_Q <- function(Y, d, rho) {
  Z <- radial_projection(Y)
  cos_simil <- cosine(t(Z))
  Q <- (1+rho^2-2*rho*cos_simil)^(-d)
  diag(Q) <- 0
  Qi <- sum(Q)
  Q_ij = Q/Qi
  return(Q_ij)
}

kl_divergence_grad <- function(Y, i, rho, d, P, cos_sim=NULL, Q=NULL) {
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
    cos_sim <- cosine(t(Z))
  if(is.null(Q))
    Q <- low_dimension_Q(Z, d, rho)
  n_minus_i <- (1:nrow(Y))[-i]
  (4*d*rho*colSums(t(
    sapply(n_minus_i, function(j) {
      return(Z[j,]/(1+rho^2-2*rho*cos_sim[i,j])*(Q[i,j]-P[i,j]))
    }))))
}


psc_sne <- function(X, d, rho_psc_list=NULL, rho=0.5, perplexity=30, num_iteration=200, 
                    initial_momentum=0.5, final_momentum=0.8, eta=200, 
                    early_exageration=4.0, colors=NULL, visualize_prog=FALSE) {
  if(d<1)
    stop("Error, d value must be greater or equal than 1")
  
  n <- nrow(X)
  p <- ncol(X)-1
  P_cond <- NULL
  
  if(is.null(rho_psc_list)) {
    res_opt <- rho_optim_bst(X, perplexity)
    rho_psc_list <- res_opt$rho_values
    P_cond <- res_opt$P
  } else {
    cosine_sim_polysphere <- cosine_polysph(X)
    P_cond <- high_dimension(x=X, rho_list=rho_psc_list, cos_sim_pol = cosine_sim_polysphere)
  }
    
  P <- symmetric_probs(P_cond)
  P = P * early_exageration
  
  total_iterations <- num_iteration + 2
  Y <- array(NA, c(n, d+1, total_iterations))
  Y[,,1] <- Y[,,2] <- r_unif_sphere(n, d+1)
  
  momentum = initial_momentum
  max_tol = 1e-3
  
  range_iterations <- seq_len(num_iteration) + 2
  for(i in range_iterations) {
    print(sprintf("Iteration %d", i-2))
    if(i>=250)
      momentum = final_momentum
    
    grad = t(simplify2array(mclapply(mc.cores = detectCores()-1, 1:n, 
                      kl_divergence_grad, Y=Y[,,i-1], rho=rho, d=d, P=P, 
                      cos_sim = cosine(t(Y[,,i-1])), 
                      Q=low_dimension_Q(Y[,,i-1], d, rho))))
    
    Y_i <- Y[,,i-1] + (eta*-grad) + momentum*(Y[,,i-1]-Y[,,i-2])
    Y[,,i] <- radial_projection(Y_i)
    
    if(i == 102) {
      P = P / early_exageration
    }
    
    if(visualize_prog && (i == 3 || (i-2) %% 5 == 0 || i == num_iteration+2)) {
      Q <- low_dimension_Q(Y[,,i], d, rho)
      C = sum(P * log(P/Q), na.rm=TRUE)
      print(sprintf("Iteration %d: objective function value is %f", i-2, C))
      visualize_iter_sol(Y, i, d, colors)
    }
  }
  return(Y)
}

visualize_iter_sol <- function(Y, i, d, colors) {
  if(is.null(colors))
    colors <- rep(1, n)
  if(d==1) {
    Y_rad <- DirStats::to_rad(Y[,,i])
    r <- 1
    theta <<- Y_rad
    plot(r*sin(theta),
         r*cos(theta),
         col=colors,
         xlim=c(-max(r),max(r)),
         ylim=c(-max(r),max(r)))
    
    polygon(max(r)*sin(seq(0,2*pi,length.out=100)),max(r)*cos(seq(0,2*pi,length.out=100)))
  } else if(d==2) {
    scatterplot3d::scatterplot3d(Y[,,i],
                                 xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
                                 color = colors)
  }
}

