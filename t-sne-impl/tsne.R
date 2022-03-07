library(stats4)
library(mvtnorm)
source("utils.R")

.high_dimension_pij <- function(X, tolerance=1e-5, perplexity=30) {
  n <- nrow(X)
  D <- x_diff(X)
  P <- matrix(0, nrow = 150, ncol = 150)
  beta <- rep(1, n)
  log_perp <- log(perplexity)
  for (i in seq_len(n)) {
    column_index <- index_except_i(i, n)
    D_i <- D[i, column_index]
    res <- entropy_beta(D_i, beta[i])
    h_star <- res$entropy
    prob_star <- res$probs
    h_diff <- h_star - log_perp
    res_opt <- binary_search_optimization(
      D_i, i, beta, h_star, prob_star,
      log_perp, tolerance
    )
    prob_star <- res_opt$probs
    beta <- res_opt$beta
    P[i, column_index] <- prob_star
  }
  print("Values of sigma for each x_i")
  print(sqrt(1 / beta))
  print(sprintf("Mean value of sigma: %01.2f", mean(sqrt(1 / beta))))
  return(P)
}

.low_dimension_qij <- function(X, P, q, num_iteration=1500, initial_momentum=0.5,
                               final_momentum=0.8, eta=500) {
  dY = replicate(q, rep(0, n))
  iY = replicate(q, rep(0, n))
  total_iterations <- num_iteration + 2
  Y <- array(NA, c(n,q,total_iterations))
  Y[,,1] <- Y[,,2] <- mvtnorm::rmvnorm(n = n, 
                                       mean = rep(0, q), 
                                       sigma = (1e-4* diag(1, q)))
  
  range_iterations <- seq_len(num_iteration) + 2
  for(i in range_iterations) {
    sum_Q <- apply(Y[,,i-1]^2, MARGIN = 1, FUN = sum)
    sum_Q <- replicate(n, sum_Q)
    num <- -2 * (Y[,,i-1] %*% t(Y[,,i-1]))
    num <- 1 / (1 + (t(num + sum_Q) + sum_Q))
    diag(num) <- 0
    Q <- num / sum(num)
    Q <- replace(Q, Q < 1e-12, 1e-12)
    
    PQ = P-Q
    
    for(j in seq_len(n)) {
      gradient = replicate(q, PQ[, j] * num[, j]) * 
        (t(replicate(n, Y[j, ,q])) - Y[,,q])
      dY[j, ] = apply(gradient, FUN=sum, MARGIN=2)
    }
    
    momentum <- initial_momentum
    if(i < 250) {
      momentum = final_momentum
    } 
    
    iY = (eta * dY) + (momentum * (Y[,,i-1]-Y[,,i-2]))
    Y[,,i] = Y[,,i-1] + iY
    
    if(i %% 10 == 0) {
      C = sum(P * log(P/Q))
      print(sprintf("Iteration %d: error is %f", i, C))
    }
    
    if(i == 100) {
      P = P / 4
    }
  }
  return(Y[,,total_iterations])
}

.execute_reg_data <- function(X) {
  high_dim_probs <- .high_dimension_pij(x, perplexity)
  high_dim_probs = .symmetric_probs(high_dim_probs)
  high_dim_probs = high_dim_probs * 4
  high_dim_probs = replace(high_dim_probs, high_dim_probs < 1e-12, 1e-12)
  
  q_matrix <- .low_dimension_qij(x, q, num_iteration, initial_momentum, 
                                 final_momentum, learning_rate)
}

execute.default <- .execute_reg_data

execute <- function(x) {
  UseMethod("execute", x)
}

tsne <- setRefClass("tsne_gen",
  fields = list(
    perplexity = "numeric", num_iteration = "numeric",
    learning_rate = "numeric", initial_momentum = "numeric",
    final_momentum = "numeric", exageration = "logical", q = "integer"
  ),
  methods = list(execute = execute)
)

print.tsne_gen <- function(obj) {
  cat(obj$perplexity, "\n")
  cat(obj$num_iteration, "\n")
  cat(obj$learning_rate, "\n")
  cat(obj$momentum, "\n")
}
