
#' @title Kullback-Leibler divergence gradient
#'
#' @description Calculates analytically the gradient of the Kullback-Leibler divergence function.
#'
#' @inheritParams low_dimension_Q
#' @inheritParams d_sph_cauchy
#' @param d target dimension to reduced the data.
#' @param P high-dimensional polyspherical Cauchy probabilities.
#' @param cos_sim an array of size \code{c(n, n, r)} with the cosine similarities in high-dimension for the polysphere \eqn{(\mathcal{S}^p)^r}. Optional parameter, default value set to \code{NULL}.
#' @param Q low-dimension spherical Cauchy probabilities. Optional parameter, default value set to \code{NULL}.
#' @return Resulting reduced data for the \eqn{i}-th observation onto the sphere \eqn{\mathcal{S}^d}
kl_divergence_grad <- function(Y, i, rho, d, P, cos_sim = NULL, Q = NULL) {
  if (i < 1 || i > nrow(Y)) {
    stop("i value not allowed. Values > 0 and <= nrow(Y)")
  }
  if (d < 1) {
    stop("d value not allowed, d >= 1")
  }
  if (rho < 0 || rho >= 1) {
    stop("rho value not allowed, d in [0, 1)")
  }
  if (nrow(Y) != nrow(P)) {
    stop("P size does not match with the num of observations")
  }
  if (d + 1 != ncol(Y)) {
    stop("d does not match with ncol(Y) - 1")
  }
  # Calculate the radial projection in case it is not in the sphere or radius 1
  Z <- radial_projection(Y)
  # Calculate the cosine similarities in case it is not passing as a parameter
  if (is.null(cos_sim)) {
    cos_sim <- reconstruct_cos_sim_mat(
      cos_sim_vec = sphunif::Psi_mat(array(Z, dim = c(nrow(Z), ncol(Z), 1)), scalar_prod = TRUE),
      n = nrow(Z)
    )
  }
  # Calculate the low dimension probabilities based on the data Z
  if (is.null(Q)) {
    Q <- low_dimension_Q(Z, rho)
  }
  # Applying the formula of the gradient: 4dp\sum_{j=1}^n[y_i'/(1+rho^2-2rho*y_i'*y_j)*(q_{ij}-p_{ij})]
  n_minus_i <- (1:nrow(Y))[-i]
  (4 * d * rho * colSums(t(
    sapply(n_minus_i, function(j) {
      return(Z[j, ] / (1 + rho^2 - 2 * rho * cos_sim[i, j]) * (Q[i, j] - P[i, j]))
    })
  )))
}

#' @title Polyspherical Cauchy SNE
#'
#' @description Calculates the polyspherical-Cauchy SNE given a data onto the poly-sphere and the reduced dimension.
#'
#' @param X an array of size \code{c(n, d + 1, r)} with the polyspherical data, where \code{n} is the number of observations, \code{d} is the dimension of each sphere, and \code{r} is the number of spheres.
#' @param d the target dimension to use for reduce the dimension of the data \code{X}.
#' @param rho_psc_list rho parameters of the high-dimensional polyspherical Cauchy probabilities (optional, default \code{NULL}).
#' @param rho rho parameter of the low-dimensional spherical Cauchy probabilities (optional, default value \code{0.5}).
#' @param perplexity parameter which says what is more important: local or global aspects (optional, default \code{30}).
#' @param num_iteration maximum number of iterations (optional, default value \code{200}).
#' @param initial_momentum first value of the momentum of the first 250 iterations.
#' @param final_momentum momentum to take into account after the 250 iteration.
#' @param eta is the learning rate of the optimization algorithm (optional, default value \code{200}).
#' @param early_exaggeration the first 100 iterations results are exaggerated, k times (optional, default value \code{4.0}).
#' @param colors list with as many elements as observations are, only valid when visualization is true (optional, default value \code{NULL}).
#' @param visualize_prog defines whether the progression plots are shown or not (optional, default value \code{FALSE}).
#' @param tol is the tolerance, when is below this value it is considered that a good solution has been obtained (optional, default value \code{1e-9}).
#' @param check whether to check or not the tolerance.
#' @param parallel_cores number of cores to use concurrently for the calculation of the gradient.
#' @return Resulting data reduced to \eqn{\mathcal{S}^d} after applying the algorithm for the total number of iterations selected.
#' @export
#' @examples
#' X <- sphunif::r_unif_sph(40, 3, 3)
#' psc_sne(X, d = 1, parallel_cores = 2)
#' psc_sne(X, d = 2, parallel_cores = 2)
psc_sne <- function(X, d, rho_psc_list = NULL, rho = 0.5, perplexity = 30, num_iteration = 200,
                    initial_momentum = 0.5, final_momentum = 0.8, eta = 200,
                    early_exaggeration = 4.0, colors = NULL, visualize_prog = FALSE,
                    tol = 1e-9, check = TRUE, parallel_cores = parallel::detectCores() - 1) {
  if (d < 1) {
    stop("Error, d value must be greater or equal than 1")
  }

  # Initializing the vectors for the objective value, errors and gradient norm
  obj_func_iter <- numeric(length = nrow(X))
  absolute_errors <- numeric(length = nrow(X))
  relative_errors <- numeric(length = nrow(X))
  gradient_norms <- numeric(length = nrow(X))

  # Sample size
  n <- nrow(X)
  # Number of dimensions
  p <- ncol(X) - 1

  # Initialize to null
  # Conditional probabilities of high-dimension
  P_cond <- NULL
  # Current data onto the reduced sphere
  Y_i <- NULL

  # Based on the rho values parameter, do different things
  if (is.null(rho_psc_list)) {
    # Calculating the rho optimal values by means of the 'rho_optim_bst' method
    res_opt <- rho_optim_bst(X, perplexity, clusterFactory(parallel_cores))
    P_cond <- res_opt$P
    rho_psc_list <- res_opt$rho_values
  } else if (is.list(rho_psc_list)) {
    # Obtaining the values from the object of the parameter
    P_cond <- rho_psc_list$P
    rho_psc_list <- rho_psc_list$rho_values
  } else {
    # Calculating the probabilities based on the rho values
    cos_sim_psh <- cosine_polysph(X)
    P_cond <- high_dimension(x = X, rho_list = rho_psc_list, cos_sim_psh = cos_sim_psh)
  }

  # Generating the high-dimension symmetric P matrix, (P_j|i + P_i|j) / 2n
  P <- symmetric_probs(P_cond)

  # Early exaggeration in the high-dimensional probabilities
  P <- P * early_exaggeration

  # Matrices to store the Y's and the Q's in each iteration
  total_iterations <- num_iteration + 2
  Y <- array(NA, c(n, d + 1, total_iterations))
  # Generate points evenly spaced
  Y[, , 1] <- Y[, , 2] <- gen_opt_sphere(n, d)
  Qs <- array(NA, c(n, n, total_iterations))
  # Generate low-dimension probabilities for the data generated
  Qs[, , 1] <- Qs[, , 2] <- low_dimension_Q(Y[, , 2], rho)

  # Initial momentum
  momentum <- initial_momentum

  # Visualizing the plots in a 4x4 grid
  old_par <- par(mfrow = c(4, 4))

  # Interval from 2 to number of iterations + 2
  range_iterations <- seq_len(num_iteration) + 2

  for (i in range_iterations) {
    # applying final momentum
    if (i >= 250) {
      momentum <- final_momentum
    }

    # Calculate the cosine similarities for the current Y solution
    Y_cos_sim <- reconstruct_cos_sim_mat(
      sphunif::Psi_mat(array(Y[, , i - 1], dim = c(nrow(Y), ncol(Y), 1)), scalar_prod = TRUE),
      nrow(Y)
    )

    # gradient of the objective function for all the observations
    grad <- t(simplify2array(parallel::mclapply(
      mc.cores = parallel_cores, 1:n,
      kl_divergence_grad, Y = Y[, , i - 1], rho = rho, d = d, P = P,
      cos_sim = Y_cos_sim,
      Q = Qs[, , i - 1]
    )))

    # gradient descent
    Y_i <- Y[, , i - 1] + (eta * -grad) + momentum * (Y[, , i - 1] - Y[, , i - 2])

    # Projecting iteration solution onto the sphere/circumference of radio 1
    Y_i <- radial_projection(Y_i)

    # Store the iteration's Y in the 3d Y's matrix
    Y[, , i] <- Y_i

    # Generate the Q matrix with the low-dimension probabilities
    Qs[, , i] <- low_dimension_Q(Y[, , i], rho)

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

    # Reverse the early exaggeration made at the beginning
    if (i == 102) {
      P <- P / early_exaggeration
    }

    # Relative error less than 1%, then break the loop
    if (check && i - 2 > 100) {
      if (relative_errors[i - 2] < tol) break
    }
  }
  # Undo the par configuration (grid 4x4)
  par(old_par)
  # Visualize the solution
  visualize_iter_sol(Y, i, d, colors)
  return(Y_i)
}

#' @title Visualize the reduction status
#'
#' @description Visualize the iteration solution in a plot for the \eqn{i}-th iteration.
#'
#' @inheritParams low_dimension_Q
#' @param i the \eqn{i}-th iteration.
#' @param d the dimension to reduce the original data.
#' @param colors defines the group colors in the plot. Optional parameter, default value set to \code{NULL}.
visualize_iter_sol <- function(Y, i, d, colors = NULL) {
  # If colors is null, set all of them to black
  if (is.null(colors)) {
    n <- nrow(Y)
    colors <- rep(1, n)
  }
  # Plot in a circumference
  if (d == 1) {
    Y_rad <- DirStats::to_rad(Y[, , i])
    r <- 1
    theta <- Y_rad
    plot(r * sin(theta),
      r * cos(theta),
      col = colors,
      xlim = c(-max(r), max(r)),
      ylim = c(-max(r), max(r)), main = paste("Iteration", i - 2)
    )

    graphics::polygon(max(r) * sin(seq(0, 2 * pi, length.out = 100)), max(r) * cos(seq(0, 2 * pi, length.out = 100)))
  }
  # Plot in an sphere
  else if (d == 2) {
    scatterplot3d::scatterplot3d(Y[, , i],
      xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
      color = colors, main = paste("Iteration", i - 2)
    )
  }
}
