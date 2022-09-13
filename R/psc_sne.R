
#' @title Kullback--Leibler divergence gradient
#'
#' @description Calculates analytically the gradient of the Kullback--Leibler
#' divergence between low- and high-dimensional pairwise probabilities.
#'
#' @inheritParams low_dimension_Q
#' @inheritParams d_sph_cauchy
#' @param d target dimension to reduce the data.
#' @param P matrix of size \code{c(n, n)} with the high-dimensional
#' polyspherical Cauchy probabilities.
#' @param cos_sim a vector of size \code{n/2} with the cosine
#' similarities in high-dimension for the polysphere \eqn{(\mathcal{S}^p)^r}.
#' The way that the cosine similarities matrix is treated makes the calculus
#' faster since it is flat in a vector object. Optional parameter, defaults to
#' \code{NULL}.
#' @param Q matrix of size \code{c(n, n)} with the low-dimension spherical
#' Cauchy probabilities. Optional parameter, defaults to \code{NULL}.
#' @return Resulting reduced data for the \eqn{i}-th observation onto the
#' sphere \eqn{\mathcal{S}^d}.
#' @examples
#' Y <- sphunif::r_unif_sph(40, 3, 1)[ , , 1]
#' X <- sphunif::r_unif_sph(40, 3, 3)
#' P <- high_dimension(X, rep(0.5, 40), num_cores = 2)
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

  # Applying the formula of the gradient for C
  # 4 * d * rho \sum_{j = 1}^n [y_i' / (1 + rho^2 - 2 * rho * y_i' y_j) *
  #                             (q_{ij} - p_{ij})]
  grad <- 4 * d * rho *
    colSums(Z[-i, ] / (1 + rho^2 - 2 * rho * cos_sim[i, -i]) *
              (Q[i, -i] - P[i, -i]))

  # Applying the formula of the gradient for C_bar
  I <- diag(rep(1, ncol(Z)))
  return(drop(grad %*% (I - tcrossprod(Z[i, ]))))

}


#' @title Polyspherical Cauchy SNE
#'
#' @description Computes polyspherical Cauchy SNE given a data onto the
#' polysphere and a dimension to reduce to.
#'
#' @param X an array of size \code{c(n, d + 1, r)} with the polyspherical data,
#' where \code{n} is the number of observations, \code{d} is the dimension of
#' each sphere, and \code{r} is the number of spheres.
#' @param d the target dimension to use for the reduced data of \code{X}.
#' @param rho_psc_list rho parameters of the high-dimensional polyspherical
#' Cauchy probabilities. Multiple types of parameters are allowed,
#' distinguishing three scenarios. The first one when the type of the parameter
#' is a list, then it contains the vector \code{rho_values} and the matrix
#' \code{P}, the second scenario when the type is a vector, then this object
#' contains the rho values, within this function the
#' \code{\link{high_dimension}} function is called to get the matrix \code{P}.
#' The last scenario that is when this object is set to \code{NULL}, i.e., the
#' \code{\link{rho_optim_bst}} function is called to get the rho values
#' (given a fixed perplexity) and the probabilities matrix. Optional parameter,
#' defaults to \code{NULL}.
#' @param rho parameter of the low-dimensional spherical Cauchy
#' probabilities. Optional, defaults to \code{0.5}.
#' @param perplexity parameter that measures the number of neighbors to
#' use when mapping between high- and low-dimension. Defaults to \code{30}.
#' @param maxit maximum number of iterations. Defaults to \code{1e3}.
#' @param initial_momentum first value of the momentum of the first \code{250}
#' iterations. Defaults to \code{0.5}.
#' @param final_momentum momentum to take into account after the \code{250}
#' iteration. Defaults to \code{0.8}.
#' @param eta is the learning rate of the optimization algorithm. Optional
#' param, defaults to \code{200}.
#' @param early_exaggeration the first \code{100} iterations results are
#' exaggerated by this factor. Optional parameter, defaults to \code{4.0}.
#' @param colors list with as many elements as observations are, only valid
#' when visualization is true. Optional parameter, defaults to \code{NULL}.
#' @param show_prog defines the number of iterations skipped when reporting
#' the progress. Defaults to \code{100}, i.e., only multiples of \code{100}
#' are reported. If \code{FALSE}, no progress is shown at all.
#' @param show_plots show convergence plots? If \code{TRUE} (default),
#' a plot is shown: after \code{2 * show_prog} iterations and at the end
#' of the search.
#' @param tol is the tolerance, when is below this value it is considered that
#' a good solution has been obtained. Defaults to \code{1e-6}.
#' @param parallel_cores number of cores to use concurrently for the
#' calculation of the gradient. Defaults to \code{parallel::detectCores() - 1},
#' that means that uses the total number of cores of the computer except one of
#' them.
#' @param init how to initialize the scores: \code{"equispaced"} (evenly spaced
#' points on the circumference/sphere), \code{"random"} (random points
#' generated uniformly), \code{"most_promising"} (best configuration
#' obtained in \code{N} differently-initialized searches run with
#' \code{maxit / 100} iterations), or a matrix. Defaults to \code{"equispaced"}.
#' @param N number of differently-initialized searches (see above). Defaults to
#' \code{10}.
#' @return A list with the following entries:
#' \itemize{
#'   \item \code{best_Y}: best configuration of scores found.
#'   \item \code{last_Y}: last configuration of scores found.
#'   \item \code{rho_psc_list}: vector or rho's.
#'   \item \code{diagnostics}: data frame with the objective function values,
#'   absolute/relative errors, gradient norms, and moment norms.
#'   \item \code{convergence}: convergence flag.
#' }
#' @details
#' When \code{init = "most_promising"}, \code{N - 1} initializations are random
#' and one is an equispaced grid.
#' @examples
#' X <- sphunif::r_unif_sph(n = 100, p = 3, M = 3)
#' X[1:50, , 1] <- rotasym::r_vMF(n = 50, mu = c(0, 0, 1), kappa = 10)
#' X[51:100, , 1] <- rotasym::r_vMF(n = 50, mu = c(0, 0, -1), kappa = 10)
#' psc <- psc_sne(X = X, d = 1, parallel_cores = 2, eta = 50,
#'                init = "most_promising")
#' psc2 <- psc_sne(X = X, d = 1, parallel_cores = 2, eta = 50,
#'                 init = psc$last_Y)
#' @export
psc_sne <- function(X, d, rho_psc_list = NULL, rho = 0.5, perplexity = 30,
                    maxit = 1e3, initial_momentum = 0.5,
                    final_momentum = 0.8, eta = 200, early_exaggeration = 4.0,
                    colors = NULL, show_prog = 100, show_plots = TRUE,
                    tol = 1e-6, parallel_cores = parallel::detectCores() - 1,
                    init = c("equispaced", "random", "most_promising")[1],
                    N = 10) {

  # Input checks
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
      (!is.numeric(maxit) || maxit <= 0) &&
      (!is.numeric(eta) || eta <= 0) &&
      (!is.numeric(parallel_cores) || parallel_cores <= 0)) {

    stop("Error, perplexity, maxit, eta, early_exaggeration and parallel_cores must be positive number")

  }

  if (!is.numeric(tol) || tol >= 1 || tol < 0) {

    stop("Error, tol parameter must be within [0, 1)")

  }

  if (N < 1) {

    stop("N has to be larger than 1")

  }

  # Initializing the vectors for the objective value, errors and gradient norm
  obj_func_iter <- numeric(length = nrow(X))
  absolute_errors <- numeric(length = nrow(X))
  relative_errors <- numeric(length = nrow(X))
  gradient_norms <- numeric(length = nrow(X))
  moment_norms <- numeric(length = nrow(X))

  # Sample size and sphere dimension within polysphere
  n <- nrow(X)

  # Conditional high-dimensional probabilities and data on the reduced sphere
  P_cond <- Y_i <- NULL

  # If show_prog is TRUE
  if (show_prog && !is.numeric(show_prog)) {

    # Set the print trace each 100 iterations
    show_prog <- 100

  }

  # Based on the rho values parameter, do different things
  if (is.null(rho_psc_list)) {

    # Calculating the rho optimal values by means of the 'rho_optim_bst' method
    if (show_prog) {

      message("Computing rho_psc_list with rho_optim_bst()")

    }
    res_opt <- rho_optim_bst(x = X, perp_fixed = perplexity,
                             num_cores = parallel_cores)
    P_cond <- res_opt$P
    rho_psc_list <- res_opt$rho_values

  } else if (is.list(rho_psc_list)) {

    # Obtaining the values from the object of the parameter
    P_cond <- rho_psc_list$P
    rho_psc_list <- rho_psc_list$rho_values

  } else {

    # Calculating the probabilities based on the rho values
    P_cond <- high_dimension(x = X, rho_list = rho_psc_list,
                             num_cores = parallel_cores)

  }

  # Generating the high-dimension symmetric P matrix, (P_j|i + P_i|j) / (2n)
  P <- symmetric_probs(P_cond)

  # Early exaggeration in the high-dimensional probabilities
  P <- P * early_exaggeration

  # Matrices to store the two previous Y's
  Y <- array(data = NA, dim = c(n, d + 1, 2))

  # Initialization
  if (is.character(init) && init == "equispaced") {

    # Generate points evenly spaced
    Y[, , 1] <- Y[, , 2] <- grid_sphere(n, d)

  } else if (is.character(init) && init == "random") {

    # Generate random points uniformly
    Y[, , 1] <- Y[, , 2] <- sphunif::r_unif_sph(n, d + 1)[, , 1]

  } else if (is.character(init) && init == "most_promising") {

    # Maximum iterations
    maxit_start <- max(round(maxit / 10), 1)

    # Run short psc_sne()'s
    start <- list()
    i_min <- 1
    obj_min <- Inf
    for (i in seq_len(N)) {

      if (show_prog) {

        message(sprintf("BEGIN start %d / %d", i, N))

      }

      # Fit and best objective explored
      start[[i]] <- psc_sne(X = X, d = d, rho_psc_list = rho_psc_list,
                            rho = rho, perplexity = perplexity,
                            maxit = maxit_start,
                            initial_momentum = initial_momentum,
                            final_momentum = final_momentum, eta = eta,
                            early_exaggeration = early_exaggeration,
                            colors = colors, show_prog = FALSE,
                            show_plots = FALSE, tol = tol,
                            parallel_cores = parallel_cores,
                            init = ifelse(i == 1, "equispaced", "random"))
      obj_i <- tail(start[[i]]$diagnostics$obj, n = 1)

      # Better solution?
      if (obj_i < obj_min) {

        obj_min <- obj_i
        i_min <- i

      }

      if (show_prog) {

        message(sprintf(paste("END start %d / %d (best: %d / %d).",
                              "Obj: %.2e (best: %.2e)"), i, N, i_min, N,
                        obj_i, obj_min))

      }

    }

    # Best solution explored
    Y[, , 1] <- Y[, , 2] <- start[[i_min]]$last_Y

  } else if (is.matrix(init)) {

    if (!(nrow(init) == n && ncol(init) == d + 1)) {


      stop("Dimension of init is not c(n, d + 1).")

    }
    Y[, , 1] <- Y[, , 2] <- init

  } else {

    stop("init must be \"equispaced\", \"random\", \"most_promising\" or a matrix.")

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

  if (show_prog && show_plots) {

      # Visualizing the plots in a 3 x 3 grid
      old_par <- par()
      on.exit(par(old_par))
      old_par <- par(mfrow = c(3, 3), mar = c(0, 0, 2, 0))

  }

  # Interval from 3 to number of iterations + 2
  range_iterations <- seq_len(maxit) + 2

  # Gradient descent
  convergence <- FALSE
  for (i in range_iterations) {

    # Applying final momentum
    if (i - 2 >= 250) {

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
      mc.cores = parallel_cores, X = seq_len(n),
      FUN = kl_divergence_grad, Y = Y[, , 2], rho = rho, d = d, P = P,
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

    # Objective function value, absolute and relative errors, gradient norm,
    # and moment norm
    obj_func_iter[i - 2] <- sum(P * log(P / Q_i), na.rm = TRUE)
    if (i > 3) {

      absolute_errors[i - 2] <- abs(obj_func_iter[i - 3] - obj_func_iter[i - 2])
      relative_errors[i - 2] <- absolute_errors[i - 2] / obj_func_iter[i - 3]

    }
    gradient_norms[i - 2] <- sqrt(sum(grad^2))
    moment_norms[i - 2] <- sqrt(sum(moment_i^2))

    # When the objective function value just calculated is smaller than the best
    if (obj_func_iter[i - 2] < best_obj_i) {

      best_i <- i - 2
      best_obj_i <- obj_func_iter[best_i]
      best_Y_i <- Y[, , 2]

    }

    # Progress trace
    if (show_prog && is.numeric(show_prog) &&
        (((i - 2) %% show_prog == 0) || (i - 2 == 1) || i - 2 == maxit)) {

      # Progress trace
      message(sprintf(paste("It: %d (best: %d); obj: %.2e (best: %.2e);",
                            "abs: %.1e; rel: %.1e; norm: %.1e; mom: %.1e"),
                      i - 2, best_i, obj_func_iter[i - 2], best_obj_i,
                      absolute_errors[i - 2], relative_errors[i - 2],
                      gradient_norms[i - 2], sqrt(sum(moment_i^2))
      ))

    }

    # Plot the current status if the first iteration, the last one or it is
    # twice the number of lines indicated with show_prog param
    if (show_prog && show_plots &&
        (i - 2 == 1 || ((i - 2) %% (show_prog * 2) == 0) || i - 2 == maxit)) {

      show_iter_sol(Y = Y[, , 2], i = i, d = d, colors = colors)

    }

    # Reverse the early exaggeration made at the beginning
    if (i - 2 == 100) {

      P <- P / early_exaggeration

    }

    # Relative error less than tol, then break the loop
    if (i - 2 > 100) {

      if (all(c(relative_errors[i - 2], absolute_errors[i - 2],
                gradient_norms[i - 2]) < tol)) {

        # Progress trace
        message(sprintf(paste("It: %d (best: %d); obj: %.2e (best: %.2e);",
                              "abs: %.1e; rel: %.1e; norm: %.1e; mom: %.1e"),
          i - 2, best_i, obj_func_iter[i - 2], best_obj_i,
          absolute_errors[i - 2], relative_errors[i - 2], gradient_norms[i - 2],
          sqrt(sum(moment_i^2))
        ))

        # Show the last configuration
        if (show_plots) {

          show_iter_sol(Y = Y[, , 2], i = i, d = d, colors = colors)

        }
        convergence <- TRUE
        break

      }
    }
  }

  # Converged?
  if (show_prog) {

    message(ifelse(
      convergence, "CONVERGENCE!",
      "**NO** CONVERGENCE. Decrease eta? Change init? Increase maxit?"))

  }

  # Return configurations and diagnostics
  return(list("best_Y" = best_Y_i, "last_Y" = Y[, , 2],
              "rho_psc_list" = rho_psc_list,
              "diagnostics" = data.frame("obj" = obj_func_iter,
                                         "abs" = absolute_errors,
                                         "rel" = relative_errors,
                                         "grad" = gradient_norms,
                                         "mom" = moment_norms),
              "convergence" = convergence))

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
#' @keywords internal
show_iter_sol <- function(Y, i, d, colors = NULL) {

  # If colors is null, set all of them to black
  if (is.null(colors)) {

    n <- nrow(Y)
    colors <- rep(1, n)

  }

  if (d == 1) {

    # Plot on a circumference
    plot(Y[, 1], Y[, 2], col = colors, xlim = c(-1, 1), ylim = c(-1, 1),
         xlab = "", ylab = "", axes = FALSE, main = paste("Iteration", i - 2))
    th <- seq(0, 2 * pi, length.out = 100)
    graphics::polygon(x = cos(th), y = sin(th))

  } else if (d == 2) {

    # Sequence from -180 to 180 by an step of 15 in radians
    seq_rad <- seq(-pi, pi, by = pi / 30)
    # Meridian calculates as theta = 0 and phi = i
    # where i is the radians
    meridian <- do.call(rbind, lapply(seq_rad, function(i) c(0, i)))
    equator <- do.call(rbind, lapply(seq_rad, function(i) c(i, pi / 2)))
    # Plot on the sphere
    sd3 <- scatterplot3d::scatterplot3d(
      Y, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
      color = colors, main = paste("Iteration", i - 2), xlab = "", ylab = "",
      zlab = "", axis = FALSE,
      pch = c('+', '-')[ifelse(sign(Y[, 2]) == 1, 1, 2)]
    )
    sd3$points3d(DirStats::to_sph(th = meridian[, 1], ph = meridian[, 2]),
                 type = "l", lty = 3)
    sd3$points3d(DirStats::to_sph(th = equator[, 1], ph = equator[, 2]),
                 type = "l", lty = 3)

  }

}
