
#' @title Kernel mean shift clustering for directional data
#'
#' @description Performs kernel mean shift clustering on \eqn{S^d} using
#' an adapted Euler algorithm and kernel density estimator.
#'
#' @param x a matrix of size \code{c(nx, d + 1)} with the initial points for
#' the Euler algorithm.
#' @param data a matrix of size \code{c(n, d + 1)} with the sample.
#' @param h bandwidth. Chosen automatically if \code{NULL} (default).
#' @param N maximum number of iterations. Defaults to \code{500}.
#' @param eps convergence tolerance. Defaults to \code{1e-3}.
#' @param tol tolerance for equality of modes. Defaults to \code{1e-1}.
#' @param keep_paths keep the ascending paths? Defaults to \code{FALSE}.
#' @param show_prog display progress?
#' @return A list with the following entries:
#' \itemize{
#'   \item \code{end_points}: end points of the Euler algorithm. A matrix of
#'   the same size as \code{x}.
#'   \item \code{cluster}: vector giving the cluster labels.
#'   \item \code{modes}: estimated modes for each cluster (sorted).
#'   \item \code{paths}: ascension paths, if \code{keep_paths = TRUE}. A list
#'   of length \code{nx} with matrices of size \code{c(np, d + 1)}, where
#'   \code{np} is at most \code{N + 1}.
#' }
#' @examples
#' # Detection of three clusters in S^2
#' data <- rbind(
#'   rotasym::r_vMF(n = 50, mu = c(0, 0, 1), kappa = 5),
#'   rotasym::r_vMF(n = 50, mu = c(0, 0, -1), kappa = 5),
#'   rotasym::r_vMF(n = 50, mu = c(1, 0, 0), kappa = 5)
#' )
#' kms <- kms_dir(x = data, data = data, keep_paths = TRUE)
#' sd3 <- scatterplot3d::scatterplot3d(data, xlim = c(-1, 1),
#'                                     ylim = c(-1, 1), zlim = c(-1, 1),
#'                                     color = kms$cluster + 1, pch = 16,
#'                                     cex.symbol = 0.5)
#' for (i in seq_len(nrow(data))) sd3$points3d(kms$paths[[i]], type = "l",
#'                                             lty = 3)
#' sd3$points3d(kms$end_points, col = kms$cluster + 1, pch = "*", cex = 2)
#'
#' # Detection of three clusters in S^1
#' data <- rbind(
#'   rotasym::r_vMF(n = 50, mu = c(0, 1), kappa = 5),
#'   rotasym::r_vMF(n = 50, mu = c(-sqrt(2), -sqrt(2)) / 2, kappa = 5),
#'   rotasym::r_vMF(n = 50, mu = c(sqrt(2), -sqrt(2)) / 2, kappa = 5)
#' )
#' kms <- kms_dir(x = data, data = data, keep_paths = TRUE)
#' plot(data, col = kms$cluster + 1, pch = 16, xlim = c(-1.5, 1.5),
#'      ylim = c(-1.5, 1.5))
#' for (i in seq_len(nrow(data))) {
#'  l <- seq(0, 1, length.out = nrow(rbind(kms$paths[[i]])))
#'  lines(sqrt(1 + l) * kms$paths[[i]], lty = 3)
#' }
#' points(sqrt(2) * kms$end_points, col = kms$cluster + 1, pch = "*", cex = 2)
#' kms
#' @export
kms_dir <- function(x, data, h = NULL, N = 500, eps = 1e-3, tol = 1e-1,
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
      if (sqrt(sum((old_y - new_y)^2)) < eps) {

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
  tree <- hclust(dist(y))
  labels <- cutree(tree, h = tol)
  modes <- sapply(seq_len(max(labels)), function(i) {
    rotasym::spherical_mean(data = y[labels == i, ])
  })

  # Return end points
  return(list(end_points = y, cluster = labels, modes = modes, paths = paths))

}


