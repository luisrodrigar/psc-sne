library(parallel)
library(doParallel)
library(lsa)

cos_sim_ij <- function(i, j) {
  (t(i) %*% j)[1]
}

## High-dimensional space

diag_3d <- function(x, k, val) {
  diag(x[, , k]) <- val
  return(x[, , k])
}

### cosine similarity

cosine_polysph <- function(X) {
  r <- dim(X)[3]
  sapply(1:r, function(x, k) cosine(t(x[, , k])), x = X, simplify = "array")
}

### Radial projection

radial_projection <- function(X) {
  t(sapply(1:nrow(X), function(i) {
    X[i, ] / norm(X[i, ], type = "2")
  }))
}

radial_projection_ps <- function(X) {
  r <- dim(X)[3]
  sapply(1:r, function(x, k) radial_projection(x[, , k]), x = X, simplify = "array")
}

### Perplexity

low_dimension_Q <- function(Y, d, rho) {
  Z <- radial_projection(Y)
  cos_simil <- cosine(t(Z))
  Q <- (1 + rho^2 - 2 * rho * cos_simil)^(-d)
  diag(Q) <- 0
  Qi <- sum(Q)
  Q_ij <- Q / Qi
  return(Q_ij)
}

kl_divergence_grad <- function(Y, i, rho, d, P, cos_sim = NULL, Q = NULL) {
  if (i < 1 || i > nrow(Y)) {
    stop(paste(
      "Error, i value not allowed. Positive values greater tha 0 and",
      "smaller or equal than the total number of observations."
    ))
  }
  if (d < 1) {
    stop("Error, d value not allowed. Positive values greater or equal than 1.")
  }
  if (rho < 0 || rho >= 1) {
    stop(paste(
      "Error, rho value not allowed, values in between 0 and 1,",
      "the last one not included."
    ))
  }
  if (nrow(Y) != nrow(P)) {
    stop("Error, the num of observations does not match with the matrix P")
  }
  if (d + 1 != ncol(Y)) {
    stop("Error, the columns of Y does not match with the value of d")
  }
  Z <- radial_projection(Y)
  if (is.null(cos_sim)) {
    cos_sim <- cosine(t(Z))
  }
  if (is.null(Q)) {
    Q <- low_dimension_Q(Z, d, rho)
  }
  n_minus_i <- (1:nrow(Y))[-i]
  (4 * d * rho * colSums(t(
    sapply(n_minus_i, function(j) {
      return(Z[j, ] / (1 + rho^2 - 2 * rho * cos_sim[i, j]) * (Q[i, j] - P[i, j]))
    })
  )))
}


psc_sne <- function(X, d, rho_psc_list = NULL, rho = 0.5, perplexity = 30, num_iteration = 200,
                    initial_momentum = 0.5, final_momentum = 0.8, eta = 200,
                    early_exageration = 4.0, colors = NULL, visualize_prog = FALSE,
                    tol = 1e-6, check = TRUE) {
  if (d < 1) {
    stop("Error, d value must be greater or equal than 1")
  }

  # Initializing the vectors for the objective value, errors and gradient norm
  obj_func_iter <- numeric(length = nrow(X))
  absolute_errors <- numeric(length = nrow(X))
  relative_errors <- numeric(length = nrow(X))
  gradient_norms <- numeric(length = nrow(X))

  n <- nrow(X)
  p <- ncol(X) - 1
  P_cond <- NULL
  Y_i <- NULL

  # Based on the rho values parameter, calculating or processing them
  if (is.null(rho_psc_list)) {
    res_opt <- rho_optim_bst(X, perplexity)
    rho_psc_list <- res_opt$rho_values
    P_cond <- res_opt$P
  } else if (class(rho_psc_list) == "list") {
    P_cond <- rho_psc_list$`P`
    rho_psc_list <- rho_psc_list$rho_values
  } else {
    cosine_sim_polysphere <- cosine_polysph(X)
    P_cond <- high_dimension(x = X, rho_list = rho_psc_list, cos_sim_pol = cosine_sim_polysphere)
  }

  # Generating the P matrix, high-dimensional probabilities (P_j|i + P_i|j) / 2n
  P <- symmetric_probs(P_cond)

  # Early exaggeration in the high-dimensional probabilities
  P <- P * early_exageration

  # Matrices to store the Y's and the Q's in each iteration
  total_iterations <- num_iteration + 2
  Y <- array(NA, c(n, d + 1, total_iterations))
  Y[, , 1] <- Y[, , 2] <- gen_opt_sphere(n, d)
  Qs <- array(NA, c(n, n, total_iterations))
  Qs[, , 1] <- Qs[, , 2] <- low_dimension_Q(Y[, , 2], d, rho)

  # Initial momentum
  momentum <- initial_momentum

  # Visualizing the plots in a 4x4 grid
  if (visualize_prog) {
    par(mfrow = c(4, 4))
  }

  # Interval from 2 to number of iterations + 2
  range_iterations <- seq_len(num_iteration) + 2

  for (i in range_iterations) {
    # apply final momentum
    if (i >= 250) {
      momentum <- final_momentum
    }

    # gradient of the objective function
    grad <- t(simplify2array(mclapply(
      mc.cores = detectCores() - 1, 1:n,
      kl_divergence_grad, Y = Y[, , i - 1], rho = rho, d = d, P = P,
      cos_sim = cosine(t(Y[, , i - 1])),
      Q = Qs[, , i - 1]
    )))

    # gradient descent
    Y_i <- Y[, , i - 1] + (eta * -grad) + momentum * (Y[, , i - 1] - Y[, , i - 2])

    # Projecting iteration solution onto the sphere/circumference of radio 1
    Y_i <- radial_projection(Y_i)

    # Store the iteration's Y in the 3d Y's matrix
    Y[, , i] <- Y_i

    # Generate the Q matrix with the low-dimension probabilities
    Qs[, , i] <- low_dimension_Q(Y[, , i], d, rho)

    # Objective func value, absolute and relative errors and the gradient norm
    obj_func_iter[i - 2] <- sum(P * log(P / Qs[, , i]), na.rm = TRUE)
    if (i > 3) {
      absolute_errors[i - 2] <- abs(obj_func_iter[i - 3] - obj_func_iter[i - 2])
      relative_errors[i - 2] <- absolute_errors[i - 2] / obj_func_iter[i - 3]
    }
    gradient_norms[i - 2] <- norm(grad, "2")

    print(sprintf(
      "Iter %d, obj %f, abs %f, rel %f, norm %f", i - 2,
      obj_func_iter[i - 2], absolute_errors[i - 2],
      relative_errors[i - 2], gradient_norms[i - 2]
    ))

    if (visualize_prog && (i == 3 || (i - 2) %% 25 == 0 || i == num_iteration + 2)) {
      visualize_iter_sol(Y, i, d, colors)
    }

    # Reverse the early exageration made at the begining
    if (i == 102) {
      P <- P / early_exageration
    }

    # Relative error less than 1%, then break the loop
    if (check && i - 2 > 1) {
      if (relative_errors[i - 2] < tol) break
    }
  }
  if (visualize_prog) {
    par(mfrow = c(1, 1))
  }
  return(Y_i)
}

gen_opt_sphere <- function(n, d) {
  Y <- NULL
  if (n < 1) {
    stop("Parameter n not valid, it has to be positive integer")
  }
  if (d == 1) {
    Y <- DirStats::to_cir(seq(0, 2 * pi, l = n + 1)[-(n + 1)]) # 0 == 2 * pi, so we exclude it
  } else if (d == 2) {
    Y <- X_Fib(n)
  } else {
    stop("Parameter d not valid, only 1 or 2")
  }
  return(Y)
}

X_Fib <- function(N) {
  phi_inv <- (1 - sqrt(5)) / 2
  N <- N / 2
  i <- (-N:N)[-N]
  lat <- asin(2 * i / (2 * N + 1))
  lon <- (2 * pi * i * phi_inv) %% (2 * pi)
  X <- cbind(cbind(cos(lon), sin(lon)) * cos(lat), sin(lat))
  return(X)
}

visualize_iter_sol <- function(Y, i, d, colors) {
  if (is.null(colors)) {
    colors <- rep(1, n)
  }
  if (d == 1) {
    Y_rad <- DirStats::to_rad(Y[, , i])
    r <- 1
    theta <<- Y_rad
    plot(r * sin(theta),
      r * cos(theta),
      col = colors,
      xlim = c(-max(r), max(r)),
      ylim = c(-max(r), max(r)), main = paste("Iteration", i - 2)
    )

    polygon(max(r) * sin(seq(0, 2 * pi, length.out = 100)), max(r) * cos(seq(0, 2 * pi, length.out = 100)))
  } else if (d == 2) {
    scatterplot3d::scatterplot3d(Y[, , i],
      xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
      color = colors, main = paste("Iteration", i - 2)
    )
  }
}
