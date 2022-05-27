library(dplyr)
library(rgl)
library(mvtnorm)
library(viridis)
library(purrr)
library(stats)

.high_dimension_pij <- function(X, tolerance=1e-5, perplexity=30) {
  n <- nrow(X)
  D <- data.matrix(stats::dist(X))
  P <- matrix(0, nrow = 150, ncol = 150)
  beta <- rep(1, n)
  entropy <- log(perplexity)
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

.low_dimension_qij <- function(X, P, q, num_iteration=500, initial_momentum=0.5,
                               final_momentum=0.8, eta=100, exageration=TRUE) {
  n <- nrow(X)
  dY = replicate(q, rep(0, n))
  iY = replicate(q, rep(0, n))
  total_iterations <- num_iteration + 2
  Y <- array(NA, c(n, q, total_iterations))
  Y[,,1] <- Y[,,2] <- mvtnorm::rmvnorm(n = n, 
                                       mean = rep(0, q), 
                                       sigma = (1e-4 * diag(1, q, q)))
  
  range_iterations <- seq_len(num_iteration) + 2
  for(i in range_iterations) {
    sum_Q <- replicate(n, rowSums(Y[,,i-1]^2))
    num <- -2 * (Y[,,i-1] %*% t(Y[,,i-1]))
    num <- 1 / (1 + (t(num + sum_Q) + sum_Q))
    diag(num) <- 0
    Q <- num / sum(num)
    Q <- replace(Q, Q < 1e-12, 1e-12)
    
    PQ = P-Q
    
    for(j in seq_len(n)) {
      gradient = replicate(q, PQ[, j] * num[, j]) * 
        (t(replicate(n, Y[j,,i-1])) - Y[,,i-1])
      dY[j, ] = apply(gradient, FUN=sum, MARGIN=2)
    }
    
    momentum <- initial_momentum
    if(i < 250) {
      momentum = final_momentum
    } 
    
    iY = -(eta * dY) + (momentum * (Y[,,i-1]-Y[,,i-2]))
    Y[,,i] = Y[,,i-1] + iY
    
    if(i %% 10 == 0) {
      C = sum(P * log(P/Q))
      print(sprintf("Iteration %d: error is %f", i, C))
    }
    
    if(i == 100 && exageration) {
      P = P / 4
    }
  }
  return(Y[,,total_iterations])
}

execute_tsne <- function(obj, X, q) {
  P <- .high_dimension_pij(X, obj@perplexity)
  P = symmetric_probs(P)
  if(obj@exageration) {
    P = P * 4
  }
  P = replace(P, P < 1e-12, 1e-12)
  
  Y <- .low_dimension_qij(X, P, q, obj@num_iteration, 
                          obj@initial_momentum, obj@final_momentum, 
                          obj@learning_rate, obj@exageration)
  return(Y)
}

setClass("tsne_gen",
         slots = list(perplexity = "numeric", num_iteration = "numeric",
                      learning_rate = "numeric", initial_momentum = "numeric",
                      final_momentum = "numeric", exageration = "logical"))

print.tsne_gen <- function(obj) {
  cat("The perplexity:", obj@perplexity, "\n")
  cat("The number of iterations:", obj@num_iteration, "\n")
  cat("The learning rate:", obj@learning_rate, "\n")
  cat("The initial momentum:", obj@initial_momentum, "\n")
  cat("The final momentum:", obj@final_momentum, "\n")
  cat("Is exageration enabled?", obj@exageration, "\n")
}

# X <- iris %>% dplyr::select(-Species) %>% as.matrix()
# obj <- new("tsne_gen", perplexity=30, num_iteration=750, learning_rate=100, 
#           initial_momentum=0.5, final_momentum=0.8, exageration=TRUE)

# res <- execute_tsne(obj, X, 2)
# plot(x=res[,1], y=res[,2], col=iris$Species)

# res <- execute_tsne(obj, X, 3)
# colors_3 <- viridis(3)
# species_colors <- ifelse(iris$Species == "setosa", colors_3[1], 
#                 ifelse(iris$Species == "versicolor", colors_3[2], colors_3[3]))
# plot3d(res[,1], res[,2], res[,3], pch = 30, col=species_colors)
# legend3d("topright", legend = c("Setosa", "Versicolor", "Virginica"), col = colors_3, pch=19)
