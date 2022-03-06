library(stats4)
source("utils.R")

.high_dimension_pij <- function(X, perplexity) {
  n <- nrow(X)
  p <- ncol(X)
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

.low_dimension_qij <- function(x) {
  return(0)
}

.execute_reg_data <- function(x) {
  p_matrix <- .high_dimension_pij(x, perplexity)
  q_matrix <- .low_dimension_qij(x)
}

execute.default <- .execute_reg_data

execute <- function(x) {
  UseMethod("execute", x)
}

tsne <- setRefClass("tsne_gen",
  fields = list(
    perplexity = "numeric", num_iteration = "numeric",
    learning_rate = "numeric", momentum = "numeric"
  ),
  methods = list(execute = execute)
)

print.tsne_gen <- function(obj) {
  cat(obj$perplexity, "\n")
  cat(obj$num_iteration, "\n")
  cat(obj$learning_rate, "\n")
  cat(obj$momentum, "\n")
}
