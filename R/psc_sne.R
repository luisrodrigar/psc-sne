
#' @title Kullback-Leibler divergence gradient
#'
#' @description Calculates analytically the gradient of the Kullback-Leibler
#' divergence between low- and high-dimensional pairwise probabilities.
#'
#' @inheritParams low_dimension_Q
#' @inheritParams d_sph_cauchy
#' @param d target dimension to reduce the data.
#' @param P matrix of size \code{c(n, n)} with the high-dimensional
#' polyspherical Cauchy probabilities.
#' @param cos_sim an array of size \code{c(n, n, r)} with the cosine
#' similarities in high-dimension for the polysphere \eqn{(\mathcal{S}^p)^r}.
#' Optional parameter, default value set to \code{NULL}.
#' @param Q matrix of size \code{c(n, n)} with the low-dimension spherical
#' Cauchy probabilities. Optional parameter, default value set to \code{NULL}.
#' @return Resulting reduced data for the \eqn{i}-th observation onto the
#' sphere \eqn{\mathcal{S}^d}
#' @examples
#' Y <- sphunif::r_unif_sph(40, 3, 1)[ , , 1]
#' X <- sphunif::r_unif_sph(40, 3, 3)
#' P <- high_dimension(X, rep(0.5, 40))
#' kl_divergence_grad(Y, 3, 0.5, 2, P)
#' cos_sim <- vec2matrix(
#'     vec = drop(sphunif::Psi_mat(array(X, dim = c(nrow(X), ncol(X), 1)),
#'                                    scalar_prod = TRUE)),
#'     n = nrow(X),
#'     diag_value = 1)
#' Q <- low_dimension_Q(Y, 0.5)
#' kl_divergence_grad(Y, 3, 0.5, 2, P, cos_sim, Q)
#' @export
kl_divergence_grad <- function(Y, i, rho, d, P, cos_sim = NULL, Q = NULL) {
  # Radial project and compute cosine similarities
  Z <- radial_projection(Y)
  if (is.null(cos_sim)) {

    cos_sim <- vec2matrix(
      vec = drop(sphunif::Psi_mat(array(Z, dim = c(nrow(Z), ncol(Z), 1)),
                                     scalar_prod = TRUE)),
      n = nrow(Z),
      diag_value = 1
    )

  }

  # Calculate the low dimension probabilities based on the data Z
  if (is.null(Q)) {

    Q <- low_dimension_Q(Z, rho)

  }

  # Applying the formula of the gradient:
  # 4 d p \sum_{j=1}^n [y_i' / (1+rho^2-2rho*y_i'*y_j) * (q_{ij}-p_{ij})]
  return(4 * d * rho *
    colSums(Z[-i, ] / (1 + rho^2 - 2 * rho * cos_sim[i, -i]) *
              (Q[i, -i] - P[i, -i])))

}


#' @title Polyspherical Cauchy SNE
#'
#' @description Calculates the polyspherical-Cauchy SNE given a data onto the
#' polysphere and the reduced dimension.
#'
#' @param X an array of size \code{c(n, d + 1, r)} with the polyspherical data,
#' where \code{n} is the number of observations, \code{d} is the dimension of
#' each sphere, and \code{r} is the number of spheres.
#' @param d the target dimension to use for reduce the dimension of the data
#' \code{X}.
#' @param rho_psc_list rho parameters of the high-dimensional polyspherical
#' Cauchy probabilities (optional, default \code{NULL}).
#' @param rho rho parameter of the low-dimensional spherical Cauchy
#' probabilities (optional, default value \code{0.5}).
#' @param perplexity parameter which says what is more important: local or
#' global aspects (optional, default \code{30}).
#' @param num_iteration maximum number of iterations (optional, default value
#' \code{200}).
#' @param initial_momentum first value of the momentum of the first \code{250}
#' iterations.
#' @param final_momentum momentum to take into account after the \code{250}
#' iteration.
#' @param eta is the learning rate of the optimization algorithm (optional,
#' default value \code{200}).
#' @param early_exaggeration the first \code{100} iterations results are
#' exaggerated by this factor (optional, default value \code{4.0}).
#' @param colors list with as many elements as observations are, only valid
#' when visualization is true (optional, default value \code{NULL}).
#' @param show_prog defines whether the progression plots are shown or
#' not (optional, default value \code{FALSE}).
#' @param tol is the tolerance, when is below this value it is considered that
#' a good solution has been obtained (optional, default value \code{1e-9}).
#' @param check whether to check or not the tolerance.
#' @param parallel_cores number of cores to use concurrently for the
#' calculation of the gradient.
#' @param init is a parameter to indicate how to proceed with the initialization
#' of the resultant reduced dimension object. There are two possible ways:
#' \code{equispaced} (evenly spaced points in the circumference/sphere) or
#' \code{random} (random points generated uniformly).
#' @return Resulting data reduced to \eqn{\mathcal{S}^d} after applying the
#' algorithm for the total number of iterations selected.
#' @examples
#' X <- sphunif::r_unif_sph(40, 3, 3)
#' psc_sne(X, d = 1, parallel_cores = 2, num_iteration = 1e4)
#' psc_sne(X, d = 2, parallel_cores = 2)
#' @export
psc_sne <- function(X, d, rho_psc_list = NULL, rho = 0.5, perplexity = 30,
                    num_iteration = 200, initial_momentum = 0.5,
                    final_momentum = 0.8, eta = 200, early_exaggeration = 4.0,
                    colors = NULL, show_prog = FALSE, tol = 1e-9,
                    check = TRUE, parallel_cores = parallel::detectCores() - 1,
                    init = c("equispaced", "random")[1]) {

  # Checks
  if (d < 1) {

    stop("Error, d value must be greater or equal than 1")

  }

  if (ncol(X) - 1 < 1) {

    stop("Error, dimension of X must be equal or greater than 2")

  }

  if (rho < 0 || rho >= 1) {

    stop("Error, high-dimensional concentration parameter, rho, must be within [0, 1)")

  }

  if (!is.null(rho_psc_list) && rlang::is_scalar_atomic(rho_psc_list)) {

    stop("Error, low-dimensional concentration parameters, rho list, must a vector")

  }

  if (!is.null(rho_psc_list) && !is.list(rho_psc_list) &&
      any(rho_psc_list >= 1 | rho_psc_list < 0)) {

    stop("Error, low-dimensional concentration parameters, rho list, must a vector")

  }

  if (!is.numeric(initial_momentum) ||
      (initial_momentum <= 0 && initial_momentum >= 1)) {

    stop("Error, initial momentum parameter must be within (0, 1)")

  }

  if (!is.numeric(final_momentum) ||
      (final_momentum <= 0 && final_momentum >= 1)) {

    stop("Error, final momentum parameter must be within (0, 1)")

  }

  if ((!is.numeric(perplexity) || perplexity <= 0) &&
      (!is.numeric(num_iteration) || num_iteration <= 0) &&
      (!is.numeric(eta) || eta <= 0) &&
      (!is.numeric(parallel_cores) || parallel_cores <= 0)) {

    stop("Error, perplexity, num_iteration, eta, early_exaggeration and parallel_cores must be positive number")

  }

  if(!is.numeric(tol) || tol >= 1 || tol < 0) {

    stop("Error, tol parameter must be within [0, 1)")

  }

  # Initializing the vectors for the objective value, errors and gradient norm
  obj_func_iter <- numeric(length = nrow(X))
  absolute_errors <- numeric(length = nrow(X))
  relative_errors <- numeric(length = nrow(X))
  gradient_norms <- numeric(length = nrow(X))

  # Sample size and sphere dimension within polysphere
  n <- nrow(X)
  p <- ncol(X) - 1

  # Conditional high-dimensional probabilities and data on the reduced sphere
  P_cond <- Y_i <- NULL

  # Based on the rho values parameter, do different things
  if (is.null(rho_psc_list)) {

    # Calculating the rho optimal values by means of the 'rho_optim_bst' method
    res_opt <- rho_optim_bst(X, perplexity, parallel_cores)
    P_cond <- res_opt$P
    rho_psc_list <- res_opt$rho_values

  } else if (is.list(rho_psc_list)) {

    # Obtaining the values from the object of the parameter
    P_cond <- rho_psc_list$P
    rho_psc_list <- rho_psc_list$rho_values

  } else {

    # Calculating the probabilities based on the rho values
    P_cond <- high_dimension(x = X, rho_list = rho_psc_list)

  }

  # Generating the high-dimension symmetric P matrix, (P_j|i + P_i|j) / (2n)
  P <- symmetric_probs(P_cond)

  # Early exaggeration in the high-dimensional probabilities
  P <- P * early_exaggeration

  # Matrices to store the two previous Y's
  Y <- array(NA, c(n, d + 1, 2))

  if(init == "equispaced") {

    # Generate points evenly spaced
    Y[, , 1] <- Y[, , 2] <- gen_opt_sphere(n, d)

  } else {

    # Generate random points uniformly
    Y[, , 1] <- Y[, , 2] <- sphunif::r_unif_sph(n, d + 1)[ , , 1]

  }

  # Initialize those object to store the best configuration
  # with the lowest value of the objective function
  best_Y_i <- NULL
  best_obj_i <- +Inf
  best_i <- NULL

  # Generate low-dimension probabilities for the data generated
  Q_i <- low_dimension_Q(Y[, , 2], rho)

  # Initial momentum
  momentum <- initial_momentum

  old_par <- NULL
  if (show_prog) {

    if (d == 2) {

      # Visualizing the plots in a 2 x 2 grid for d = 2
      old_par <- par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))

    } else if (d == 1) {

      # Visualizing the plots in a 4 x 4 grid for d = 1
      old_par <- par(mfrow = c(4, 4), mar = c(1.5, 1.5, 1.5, 1.5))

    }

  }

  # Interval from 2 to number of iterations + 2
  range_iterations <- seq_len(num_iteration) + 2

  for (i in range_iterations) {

    # Applying final momentum
    if (i >= 250) {

      momentum <- final_momentum

    }

    # Calculate the cosine similarities for the current Y solution
    Y_cos_sim <- vec2matrix(
      vec = drop(sphunif::Psi_mat(Y[, , 2, drop = FALSE], scalar_prod = TRUE)),
      n = nrow(Y),
      diag_value = 1
    )

    # Gradient of the objective function for all the observations
    grad <- t(simplify2array(parallel::mclapply(
      mc.cores = parallel_cores, 1:n,
      kl_divergence_grad, Y = Y[, , 2], rho = rho, d = d, P = P,
      cos_sim = Y_cos_sim, Q = Q_i
    )))

    # Moment
    moment_i <- momentum * (Y[, , 2] - Y[, , 1])

    # Gradient descent
    Y_i <- Y[, , 2] + (eta * -grad) + moment_i


    # Projecting iteration solution onto the sphere/circumference of radio 1
    Y_i <- radial_projection(Y_i)

    # Store the iteration's Y in the 3d Y's matrix
    Y[, , 1] <- Y[, , 2]
    Y[, , 2] <- Y_i

    # Generate the Q matrix with the low-dimension probabilities
    Q_i <- low_dimension_Q(Y_i, rho)

    # Objective func value, absolute and relative errors and the gradient norm
    obj_func_iter[i - 2] <- sum(P * log(P / Q_i), na.rm = TRUE)
    if (i > 3) {

      absolute_errors[i - 2] <- abs(obj_func_iter[i - 3] - obj_func_iter[i - 2])
      relative_errors[i - 2] <- absolute_errors[i - 2] / obj_func_iter[i - 3]

    }
    gradient_norms[i - 2] <- norm(grad, "2")

    # When the objective function value just calculated is smaller than the best
    if (obj_func_iter[i - 2] < best_obj_i) {

      best_i <- i - 2
      best_obj_i <- obj_func_iter[best_i]
      best_Y_i <- Y[ , , 2]

    }

    # Progress
    cat(sprintf(
      "It: %d; obj: %.3e; abs: %.3e; rel: %.3e; norm: %.3e; mom: %.3e;\nbest it: %d; best obj: %.3e\n",
      i - 2, obj_func_iter[i - 2], absolute_errors[i - 2],
      relative_errors[i - 2], gradient_norms[i - 2], norm(moment_i, "2"),
      best_i, best_obj_i
    ))

    if (show_prog && (i == 3 || (i - 2) %% 25 == 0 || i == num_iteration + 2)) {

      show_iter_sol(Y[ , , 2], i, d, colors)

    }

    # Reverse the early exaggeration made at the beginning
    if (i == 102) {

      P <- P / early_exaggeration

    }

    # Relative error less than tol, then break the loop
    if (check && i - 2 > 100) {

      if (all(c(relative_errors[i - 2], absolute_errors[i - 2],
                gradient_norms[i - 2]) < tol)) {

        break

      }

    }
  }

  if (show_prog) {

    # Visualize the solution
    par(old_par)
    show_iter_sol(best_Y_i, best_i + 2, d, colors)

  }

  return(list("best_Y" = best_Y_i, "last_Y" = Y[, , 2],
              "obj_seq" = obj_func_iter, "norm_seq" = gradient_norms))

}


#' @title Visualize the reduction status
#'
#' @description Visualize the iteration solution in a plot for the
#' \eqn{i}-th iteration.
#'
#' @inheritParams low_dimension_Q
#' @param i the \eqn{i}-th iteration.
#' @param d the dimension to reduce the original data.
#' @param colors defines the group colors in the plot. Optional parameter,
#' default value set to \code{NULL}.
show_iter_sol <- function(Y, i, d, colors = NULL) {

  # If colors is null, set all of them to black
  if (is.null(colors)) {

    n <- nrow(Y)
    colors <- rep(1, n)

  }

  # Plot on a circumference
  if (d == 1) {

    Y_rad <- DirStats::to_rad(Y)
    r <- 1
    theta <- Y_rad
    plot(r * sin(theta),
      r * cos(theta),
      col = colors,
      xlim = c(-max(r), max(r)),
      ylim = c(-max(r), max(r)), main = paste("Iteration", i - 2)
    )
    graphics::polygon(max(r) * sin(seq(0, 2 * pi, length.out = 100)),
                      max(r) * cos(seq(0, 2 * pi, length.out = 100)))

  }

  # Plot on an sphere
  else if (d == 2) {

    scatterplot3d::scatterplot3d(
      Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
      color = colors, main = paste("Iteration", i - 2)
    )

  }

}
