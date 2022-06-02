###########################################
##      spherical Cauchy Perplexity      ##
###########################################
library(parallel)
library(doParallel)
library(lsa)
library(stats)
library(optimParallel)

to_perplexity_P <- function(x, i, rho) {
  if(i < 1 || i > nrow(x))
    stop("The indexes i not valid, must be >= 1 and <= nrow(x)")
  if(!rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be an scalar")
  total_P_i <- sum(P_i_psc(x=x, i=i, rho_list=rep(rho, nrow(x))))
  Picondj <- psc_cond_given_i(x=x, i=i, rho_list=rep(rho, nrow(x)), 
                              total_P_i = total_P_i)
  entropy <- function(j) {
    return(Picondj[j] * log2(Picondj[j]))
  }
  perp <- 2^(-sum(sapply(seq_len(nrow(x))[-i], entropy)))
  Picondj <- Picondj
  return(list(perp=perp, Picondj=Picondj))
}

to_perplexity <- function(x, i, rho, cos_sim_ps=NULL) {
  if(i < 1 || i > nrow(x))
    stop("The indexes i not valid, must be >= 1 and <= nrow(x)")
  if(!rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be an scalar")
  if(!is.null(cos_sim_ps) && length(dim(cosine_polysph(x)))!=3)
    stop("Parameter cos_sim_ps must be a 3d-array")
  if(is.null(cos_sim_ps))
    cos_sim_ps <- cosine_polysph(x)
  Picondj <- high_dimension(x, rep(rho, nrow(x)), cos_sim_ps)
  entropy <- function(j) {
    (Picondj[i,j] * log2(Picondj[i,j]))
  }
  return(2^(-sum(sapply(seq_len(nrow(x))[-i], entropy))))
}

# check <- function(l) max(sapply(l, function(y) max(abs(l[[1]] - y)))) < 1e-7
# microbenchmark::microbenchmark(
#   to_perplexity_P(x, 1, 0.5),
#   to_perplexity(x, 1, 0.5),
#   check = check
# )
# 
# Unit: milliseconds
# expr       min        lq      mean    median        uq       max neval
# to_perplexity_P(x, 1, 0.5)  16.69372  17.88449  24.63309  18.65393  27.29528  234.1621   100
# to_perplexity(x, 1, 0.5) 517.03870 583.15446 625.63860 613.75799 638.68714 1649.2364   100

### matrix way

to_perp_scalar <- function(x, rho) {
  if(length(dim(x))!=3)
    stop("Dataset 'x' must be a 3d-array")
  if(rlang::is_scalar_atomic(rho))
    stop("Parameter rho must be a vector")
  return(sapply(1:nrow(x), function(i, x, rho) to_perplexity_P(x, i, rho[i])$perp, x=x, rho))
}

to_perp <- function(x, rho_list, cos_sim_ps=NULL) {
  if(length(dim(x))!=3)
    stop("Dataset 'x' must be a 3d-array")
  if(rlang::is_scalar_atomic(rho_list))
    stop("Parameter rho_list must be a vector")
  if(!is.null(cos_sim_ps) && length(dim(cos_sim_ps))!=3)
    stop("Parameter cos_sim_ps must be a 3d-array")
  if(is.null(cos_sim_ps))
    cos_sim_ps <- cosine_polysph(x)
  P <- high_dimension(x, rho_list, cos_sim_ps)
  op <- P*log2(P)
  diag(op) <- 0
  return(2^(-rowSums(op)))
}

# microbenchmark::microbenchmark(
#   to_perplexity_P(x, 1, 0.5),
#   to_perp(x, rep(.5, nrow(x)), cosine_polysph(x))
# )
# 
# Unit: milliseconds
# expr      min       lq      mean    median        uq      max neval
# to_perplexity_P(x, 1, 0.5) 16.52357  18.1450  26.02556  19.03716  32.80008 175.8253   100
# to_perp(x, rep(0.5, nrow(x)), cos_sim_sph) 96.76836 123.6451 152.28035 142.41592 149.62526 381.0817   100

# optimize rho based on the perplexity

### scalar calculus

# Time difference of 10.21872 mins

rho_optim_serial <- function(x, perplexity) {
  n <- nrow(x)
  d <- (ncol(x)-1)
  start_time <- Sys.time()
  rho_opt <- sapply(1:n, function(i){
    stats::optim(par = 0.5, 
          fn = function(rho) {
            res <- (to_perplexity_P(x, i, rho)$perp - perplexity)^2
            ifelse(is.finite(res), res, 1e6)
          },
          method="L-BFGS-B", lower=0, upper=.9999)$par
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  rho_opt <- simplify2array(rho_opt)
  return(rho_opt)
}

# Time difference of 18.12201 mins

rho_optim_parallel <- function(x, perplexity) {
  n <- nrow(x)
  d <- (ncol(x)-1)
  num_cores <- detectCores()-1
  cl <- clusterFactory(num_cores, 'log.txt')
  start_time <- Sys.time()
  rho_opt <- parLapply(cl, 1:n, function(i){
    print(i)
    stats::optim(par = 0.5, 
                 fn = function(rho) {
                   res <- (to_perplexity_P(x, i, rho)$perp - perplexity)^2
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

# > system.time(rho_optim_parallel(spokes, 22))
# Time difference of 29.04886 mins
# user   system  elapsed 
# 3.812    5.742 1743.000 

rho_optim_parall <- function(x, perplexity) {
  n <- nrow(x)
  d <- (ncol(x)-1)
  num_cores <- detectCores()-1
  cl <- clusterFactory(num_cores, 'log.txt')
  setDefaultCluster(cl=cl)
  start_time <- Sys.time()
  rho_opt <- sapply(1:n, FUN = function(i) {
    print(i)
    optimParallel(par = 0.5, 
                  fn = function(rho) {
                    print(rho)
                    res <- to_perplexity_P(x, i, rho)$perp - perplexity
                    ifelse(is.finite(res), res, 1e6)
                  },
                  lower=0, upper=.9999)$par
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  setDefaultCluster(cl=NULL)
  stopCluster(cl)
  return(rho_opt)
}

# system.time(rho_optim_parall(spokes, 22))
# Time difference of 1.281991 hours
# user   system  elapsed 
# 23.811   15.432 4615.212 

## matricidal calculus, cost function
## Time difference of 14.22428 mins
## 50 rows 25 spheres

rho_optim_ineff <- function(x, perplexity) {
  n <- nrow(x)
  start_time <- Sys.time()
  cosine_polysph <- cosine_polysph(x)
  rho_opt <- sapply(1:n, function(i){
    stats::optim(par = rep(0.5, n), 
          fn = function(rho) {
            print(rho)
            dif <- to_perp(x, rho, cosine_polysph) - perplexity
            res <- t(dif) %*% dif
            ifelse(is.finite(res), res, 1e6)
            },
          method="L-BFGS-B", lower=rep(0,n), upper=rep(.9999,n))$par
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  return(rho_opt)
}

## efficient ways to calculate the perplexity

## parallel and matrix calculus 
## Time difference of 14.56007 mins
## 50 rows 25 spheres

rho_optimize_1 <- function(x, perplexity) {
  n <- nrow(x)
  num_cores <- detectCores()-1
  cl <- makeForkCluster(num_cores, outfile="log.txt")
  start_time <- Sys.time()
  cosine_polysph <- cosine_polysph(x)
  rho_opt <- parLapply(cl, 1:n, function(i){
    stats::optim(par = 0.5, 
          fn = function(rho) {
            print(rho)
            res <- (to_perplexity(x = x, i = i, rho=rho, cosine_polysph) - perplexity)^2
            ifelse(is.finite(res), res, 1e6)
          },
          method="L-BFGS-B", lower=0, upper=.9999)$par
    print(paste("Finish observation:",i))
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  rho_opt <- simplify2array(rho_opt)
  stopCluster(cl) # Time difference of 5.836471 mins
  return(rho_opt)
}




rho_optimize_2 <- function(x, perplexity) {
  n <- nrow(x)
  num_cores <- detectCores()-1
  cl <- makeForkCluster(num_cores)
  start_time <- Sys.time()
  cosine_polysph <- cosine_polysph(x)
  rho_opt <- parLapply(cl, 1:n, function(i){
    stats::optim(par = 0.5, 
          fn = function(rho) {
            res <- (to_perplexity(x = x, i = i, rho=rho) - perplexity)^2
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

# check <- function(l) max(sapply(l, function(y) max(abs(l[[1]] - y)))) < 1e-7
# microbenchmark::microbenchmark(
#  rho_optimize(xxx, 15),
#  rho_optimize_2(xxx, 15),
#  check = check
#)

# Unit: milliseconds
# expr        min         lq       mean     median         uq        max neval
# rho_optimize(xxx, 15)   64.95212   74.77625   89.47905   84.41412   95.07279   200.1545   100
# rho_optimize_2(xxx, 15) 7582.09367 7729.17412 8035.37894 7880.25502 8186.76657 11580.5412   100

#check <- function(l) max(sapply(l, function(y) max(abs(l[[1]] - y)))) < 1e-7
#microbenchmark::microbenchmark(
#  rho_optimize(xxx, 15),
#  rho_optim_ineff(xxx, 15),
#  check = check
#)

### Binary search tree

bin_search <- function(perp_diff, rho, i, rho_min, rho_max) {
  if (perp_diff > 0 || is.na(perp_diff)) {
    rho_min <- rho
  } else {
    rho_max <- rho
  }
  rho <- (rho_min + rho_max) / 2
  return(list(rho = rho, min = rho_min, max = rho_max))
}

rho_optim_i_bst <- function(x, i, perp_fixed, tolerance = 1e-3, rho = NULL, 
                            max_tries = 20){
  if(is.null(rho))
    rho=0.5
  rho_min <- 0
  rho_max <- 0.9999
  tries <- 0
  perp_diff <- NA
  
  res <- to_perplexity_P(x, i, rho)
  perp_star <- res$perp
  prob_star <- res$Picondj
  
  if(is.na(perp_star))
    perp_diff <- -Inf
  else 
    perp_diff <- perp_star - perp_fixed
  while((is.na(perp_diff) || abs(perp_diff) > tolerance) && tries <= max_tries) {
    rho_opt <- bin_search(perp_diff, rho, i, rho_min, rho_max)
    rho <- rho_opt$rho
    rho_min <- rho_opt$min
    rho_max <- rho_opt$max
    
    res_loop <- to_perplexity_P(x, i, rho)
    perp_star <- res_loop$perp
    prob_star <- res_loop$Picondj
    if(is.na(perp_star))
      perp_diff <- -Inf
    else
      perp_diff <- perp_star - perp_fixed
    
    tries <- tries + 1
  }
  return(list(rho=rho, prob_i=prob_star))
}

rho_optim_bst <- function(x, perp_fixed){
  n <- nrow(x)
  num_cores <- detectCores()-1
  cl <- clusterFactory(num_cores)
  start_time <- Sys.time()
  res_opt <- parLapply(cl, 1:n, function(i){
    rho_optim_i_bst(x, i, perp_fixed)
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  stopCluster(cl)
  rho_values <- sapply(1:length(res_opt), function(i) res_opt[[i]]$rho)
  P <- t(sapply(1:length(res_opt), function(i) res_opt[[i]]$prob_i))
  return(list(rho_values=rho_values, P=P))
}

# > rho_optim_bst(x_array, 22)
# Time difference of 45.21138 secs
# with optim it took 2.58 min in the best case


# > aa <- rho_optim_bst(spokes, 22)
# Time difference of 5.758403 mins
# with optim func it took 28 min approx.


