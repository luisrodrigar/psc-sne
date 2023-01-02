
#' @title Kernel mean shift clustering for directional data
#'
#' @description Performs kernel mean shift clustering on \eqn{S^d} using
#' an adapted Euler algorithm and kernel density estimator.
#'
#' @param data a matrix of size \code{c(n, d + 1)} with the sample.
#' @param x a matrix of size \code{c(nx, d + 1)} with the initial points for
#' the Euler algorithm. Defaults to \code{data}.
#' @param h bandwidth. Chosen automatically if \code{NULL} (default).
#' @param N maximum number of iterations. Defaults to \code{500}.
#' @param eps convergence tolerance. Defaults to \code{1e-3}.
#' @param tol tolerance for equality of modes. Defaults to \code{1e-1}.
#' @param keep_paths keep the ascending paths? Defaults to \code{FALSE}.
#' @param show_prog display progress? Defaults to \code{TRUE}.
#' @return A list with the following entries:
#' \itemize{
#'   \item \code{end_points}: end points of the Euler algorithm. A matrix of
#'   the same size as \code{x}.
#'   \item \code{cluster}: vector giving the cluster labels.
#'   \item \code{modes}: estimated modes for each cluster (sorted).
#'   \item \code{paths}: ascension paths, if \code{keep_paths = TRUE}. A list
#'   of length \code{nx} with matrices of size \code{c(np, d + 1)}, where
#'   \code{np} is at most \code{N + 1}.
#'   \item \code{tree}: internal hierarchical clustering tree used to merge
#'   modes.
#'   \item \code{h}: used bandwidth.
#' }
#' @examples
#' # Detection of three clusters in S^2
#' samp <- rbind(
#'   rotasym::r_vMF(n = 50, mu = c(0, 0, 1), kappa = 5),
#'   rotasym::r_vMF(n = 50, mu = c(0, 0, -1), kappa = 5),
#'   rotasym::r_vMF(n = 50, mu = c(1, 0, 0), kappa = 5)
#' )
#' kms <- kms_dir(data = samp, keep_paths = TRUE)
#' sd3 <- scatterplot3d::scatterplot3d(samp, xlim = c(-1, 1),
#'                                     ylim = c(-1, 1), zlim = c(-1, 1),
#'                                     color = kms$cluster + 1, pch = 16,
#'                                     cex.symbol = 0.5)
#' for (i in seq_len(nrow(samp))) sd3$points3d(kms$paths[[i]], type = "l",
#'                                             lty = 3)
#' sd3$points3d(kms$end_points, col = kms$cluster + 1, pch = "*", cex = 2)
#'
#' # Detection of three clusters in S^1
#' samp <- rbind(
#'   rotasym::r_vMF(n = 50, mu = c(0, 1), kappa = 5),
#'   rotasym::r_vMF(n = 50, mu = c(-sqrt(2), -sqrt(2)) / 2, kappa = 5),
#'   rotasym::r_vMF(n = 50, mu = c(sqrt(2), -sqrt(2)) / 2, kappa = 5)
#' )
#' kms <- kms_dir(data = samp, keep_paths = TRUE)
#' plot(samp, col = kms$cluster + 1, pch = 16, xlim = c(-1.5, 1.5),
#'      ylim = c(-1.5, 1.5))
#' for (i in seq_len(nrow(samp))) {
#'  l <- seq(0, 1, length.out = nrow(rbind(kms$paths[[i]])))
#'  lines(sqrt(1 + l) * kms$paths[[i]], lty = 3)
#' }
#' points(sqrt(2) * kms$end_points, col = kms$cluster + 1, pch = "*", cex = 2)
#' kms
#' @export
kms_dir <- function(data, x = data, h = NULL, N = 500, eps = 1e-3, tol = 1e-1,
                    keep_paths = FALSE, show_prog = TRUE) {

  # Check dimensions
  nx <- nrow(x)
  d <- ncol(x) - 1
  n <- nrow(data)
  stopifnot(d == ncol(data) - 1)

  # Get upscaled EMI bandwidth for derivative
  if (is.null(h)) {

    h <- DirStats::bw_dir_emi(data = data)$h_opt *
      n^(1 / (d + 4)) * n^(-1 / (d + 6))

  }

  # Move ahead on kms
  step_ahead <- function(y) {

    # Projected gradient
    dkde <- numDeriv::grad(func = function(z) {
      DirStats::kde_dir(x = z / sqrt(sum(z^2)), data = data, h = h)
    }, x = y) # TODO: replace with analytical computation
    kde <- DirStats::kde_dir(x = y, data = data, h = h)
    eta <- dkde / kde

    # Advance and normalize
    z <- y + h^2 * eta
    z <- z / sqrt(sum(z^2))
    return(z)

  }

  # Presave space for the ascending paths
  if (keep_paths) {

    paths <- lapply(seq_len(nx), function(i) {
      matrix(x[i, ], nrow = N + 1, ncol = d + 1, byrow = TRUE)
      })

  } else {

    paths <- NULL

  }

  # Iterative procedure
  y <- x
  if (show_prog) pb <- txtProgressBar(style = 3)
  for (i in seq_len(nx)) {

    for (k in 1:N) {

      # Euler advance
      old_y <- y[i, ]
      new_y <- step_ahead(y = old_y)

      # Convergence?
      if (acos(sum(old_y * new_y)) < eps) {

        # Trim the saved path to the last iteration
        if (keep_paths) paths[[i]] <- paths[[i]][1:k, ]
        break

      } else {

        # Update
        y[i, ] <- new_y
        if (keep_paths) paths[[i]][k + 1, ] <- new_y

      }

    }

    # Progress?
    if (show_prog) setTxtProgressBar(pb, i / nx)

  }

  # Cluster end points
  tree <- hclust(acos(1 - 0.5 * dist(y)^2))
  labels <- cutree(tree, h = tol)
  modes <- sapply(seq_len(max(labels)), function(i) {
    rotasym::spherical_mean(data = y[labels == i, , drop = FALSE])
  })
  modes <- t(modes)

  # Return end points
  return(list(end_points = y, cluster = labels, modes = modes, paths = paths,
              tree = tree, h = h))

}

#' @title Estimation fo the eta parameter
#'
#' @description Auxiliar function to run kernel mean shift clustering with scalar values.
#' Calculate the estimation of the parameter eta.
#'
#' @param x a matrix of size \code{c(nx, d + 1)} with the initial points.
#' @param h bandwith parameter
#' @param data a matrix of size \code{c(n, d + 1)} with the sample.
#'
#' @return estimation of the eta parameter.
eta_hat <- function(x, h, data) {
  h^2 * ks::kdde(
    x = data, h = h, eval.points = x, deriv.order = 1,
    supp = 50, binned = FALSE
  )$estimate /
    ks::kde(
      x = data, h = h, eval.points = x, supp = 50,
      binned = FALSE
    )$estimate
}

#' @title Compute the euler value
#'
#' @description Auxiliar function to run kernel mean shift clustering with scalar values.
#' Calculate the estimation of the euler value.
#'
#' @param x a matrix of size \code{c(nx, d + 1)} with the initial points.
#' @param N. Defaults to \code{1e6}.
#' @param verbose. Defaults to \code{FALSE}.
#' @param epsilon. Defaults to \code{1e-5}.
#' @param eta_hat_efic eta parameter value.
#'
#' @return euler value for the given parameters.
euler <- function(x, N = 1e6, verbose = FALSE, epsilon = 1e-5, eta_hat_efic) {
  for (i in 1:N) {
    x_new <- sdetorus::toPiInt(x + eta_hat_efic(x = sdetorus::toPiInt(x)))
    if (verbose) print(x_new)
    if (abs(sdetorus::toPiInt(x_new - x)) < epsilon) {
      if (verbose) cat("Stopped at", i, "\n")
      break
    }
    x <- x_new
  }
  return(x)
}

#' @title Kernel mean shift clustering for linear data
#'
#' @description Computes the kernel mean shift clustering for scalar values.
#'
#' @param x a matrix of size \code{c(nx, d + 1)} with the initial points.
#' @param original_clusters. Vectors with the label associated to each observation.
#' Defaults to \code{NULL}.
#' @param h. Defaults to \code{NULL}. If it is not filled in, it is set
#' automatically to \code{ks::hpi(x, deriv.order = 1)}.
#'
#' @return A list with the following entries:
#' \itemize{
#'   \item \code{cluster}: vector giving the cluster labels.
#'   \item \code{unique_modes}: estimated modes for each cluster (sorted).
#'   \item \code{h}: used bandwidth.
#'   \item \code{position_antimodes}: indexes with the location of the antimodes.
#'   \item \code{labels_rle_values}: values of the label's rle
#'   \item \code{x_kms}: eval points
#' }
#' @export
kms_linear <- function(x, original_clusters = NULL, h = NULL) {
  # Hack to have periodicity
  samp <- c(x - 2 * pi, x, x + 2 * pi)

  if (is.null(h)) {
    h <- ks::hpi(x, deriv.order = 1)
  }

  # Speedup hack for kernel mean shift clustering
  x_spline <- seq(-1.35 * pi, 1.35 * pi, by = 0.005)
  y_spline <- sapply(x_spline, function(xx) {
    eta_hat(xx, data = samp, h = h)
  })
  eta_hat_spline <- splinefun(x = x_spline, y = y_spline)

  # Run kernel mean shift clustering
  x_kms <- seq(-1.25 * pi, 1.25 * pi, by = 0.001)
  kms <- numeric(length(x_kms))
  N <- length(x_kms)
  pb <- txtProgressBar(style = 3)
  for (i in 1:N) {
    kms[i] <- euler(x = x_kms[i], eta_hat_efic = eta_hat_spline)
    setTxtProgressBar(pb, i / N)
  }
  cat("\n")
  kms <- sdetorus::toPiInt(kms)

  # Obtain modes for data, and cluster them
  unique_modes <- sort(unique(round(kms, 3)), decreasing = TRUE)
  if (!is.null(original_clusters)) {
    unique_modes <- replace(unique_modes, c(2, 3), unique_modes[c(3, 2)])
  }

  arrival_modes <- sapply(x, euler, eta_hat_efic = eta_hat_spline)
  labels_modes <- match(
    x = round(arrival_modes, 3),
    table = unique_modes
  )

  # Cluster kms
  labels_kms <- match(x = round(kms, 3), table = unique_modes)

  if (!is.null(original_clusters)) {
    # Check matching with original_clusters
    tab <- table(
      "assigned" = labels_modes,
      "true" = original_clusters
    )
    print(tab)
    print(sprintf("Number of misclassified observations = %d", sum(tab) - sum(diag(tab))))
    print(sprintf("Classification rate = %f", sum(diag(tab)) / sum(tab)))
  }

  # Obtain relevant points
  labels_rle <- rle(labels_kms)
  positions_antimodes <- c(1, cumsum(labels_rle$lengths))

  # Return end points
  return(list(cluster = labels_kms, unique_modes = unique_modes,
              h = h, positions_antimodes = positions_antimodes,
              labels_rle_values = labels_rle$values, x_kms = x_kms))
}

#' @title Plot the density chart of the linear data variable
#'
#' @description Produces a density plot for the given parameters.
#'
#' @param x a matrix of size \code{c(nx, d + 1)} with the initial points.
#' @param original_clusters. Vectors with the label associated to each observation.
#' Defaults to \code{NULL}.
#' @param h. Defaults to \code{NULL}. If it is not filled in, the kernel mean shift is
#' calculated with h equal to \code{ks::hpi(x, deriv.order = 1)}.
#' @export
plot_kde <- function(x, original_clusters = NULL, kms_data = NULL, h = NULL) {
  # Hack to have periodicity
  samp <- c(x - 2 * pi, x, x + 2 * pi)

  if (is.null(kms_data)) {
    # Compute the kernel mean shift
    kms_data <- kms_linear(x, original_clusters, h)
  }

  unique_modes <- kms_data$unique_modes
  total_modes <- length(unique_modes)
  positions_antimodes <- kms_data$positions_antimodes
  labels_rle_values <- kms_data$labels_rle_values
  h <- kms_data$h
  x_kms <- kms_data$x_kms

  # Plot curve and axes
  plot(ks::kde(samp, h = h, gridsize = 1e3),
       xlim = c(-pi, pi),
       axes = FALSE, xlab = "", ylab = "",
       ylim = c(0, 0.1)
  )
  sdetorus::torusAxis(1)
  axis(2)
  abline(h = 0, lty = 3)

  # Colors for the plot
  col <- rainbow(total_modes + 2)

  # Plot modes
  for (i in seq_along(unique_modes)) {
    segments(
      x0 = unique_modes[i], y0 = 0, x1 = unique_modes[i],
      y1 = ks::kde(
        x = samp, h = h, eval.points = unique_modes[i],
        binned = FALSE
      )$estimate,
      col = col[i], lwd = 2
    )
  }

  # Plot domains of attraction
  kde <- ks::kde(x = samp, h = h, eval.points = x_kms, binned = FALSE)
  for (k in seq_along(labels_rle_values)) {
    begin <- positions_antimodes[k]
    end <- positions_antimodes[k + 1]
    polygon(
      x = c(
        kde$eval.points[begin],
        kde$eval.points[begin:end],
        kde$eval.points[end]
      ),
      y = c(0, kde$estimate[begin:end], 0),
      col = rainbow(total_modes + 2, alpha = 0.15)[labels_rle_values[k]],
      border = NA
    )
  }

  lines(kde$eval.points, kde$estimate)
  abline(v = x_kms[positions_antimodes], lty = 3)

  if (!is.null(original_clusters)) {
    # Plot rug
    for (i in unique(original_clusters)) {
      rug(samp[rep(original_clusters, 3) == i], col = col[i], ticksize = 0.03)
    }
  }
}
