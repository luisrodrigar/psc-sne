
#' @title Kernel mean shift clustering for directional data
#'
#' @description Performs kernel mean shift clustering on \eqn{S^d} using
#' an adapted Euler algorithm and kernel density estimator.
#'
#' @param data a matrix of size \code{c(n, d + 1)} with the sample.
#' @param x a matrix of size \code{c(nx, d + 1)} with the initial points for
#' the Euler algorithm. Defaults to \code{data}.
#' @param h bandwidth.
#' @param N maximum number of iterations. Defaults to \code{500}.
#' @param eps convergence tolerance. Defaults to \code{1e-3}.
#' @param tol tolerance for equality of modes. Defaults to \code{1e-1}.
#' @param keep_paths keep the ascending paths? Defaults to \code{FALSE}.
#' @param show_prog display progress? Defaults to \code{TRUE}.
#' @param init_clusters. Array with the label associated to each observation.
#' Defaults to \code{NULL}.
#' @param is_cut. Boolean value that informs wether the tree can be cut or not.
#' Defaults to \code{TRUE}.
#' @param is_numDeriv. Boolean value that choose the gradient computation method.
#' @return A list with the following entries:
#' \itemize{
#'   \item \code{end_points}: end points of the Euler algorithm. A matrix of
#'   the same size as \code{x}.
#'   \item \code{cluster}: vector giving the cluster labels.
#'   \item \code{modes}: estimated modes for each cluster (sorted).
#'   \item \code{antimodes}: estimated position of the antimodes.
#'   \item \code{paths}: ascension paths, if \code{keep_paths = TRUE}. A list
#'   of length \code{nx} with matrices of size \code{c(np, d + 1)}, where
#'   \code{np} is at most \code{N + 1}.
#'   \item \code{tree}: internal hierarchical clustering tree used to merge
#'   modes.
#'   \item \code{h}: used bandwidth.
#'   \item \code{eval.points}: initial points for the Euler algorithm.
#'   \item \code{labels_rle}: labels rle
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
kms_dir <- function(data, x = data, h, N = 500, eps = 1e-3, tol = 1e-1,
                    keep_paths = FALSE, show_prog = TRUE, init_clusters = NULL,
                    is_cut = TRUE, is_numDeriv = FALSE) {

  # Check dimensions
  nx <- nrow(x)
  d <- ncol(x) - 1
  n <- nrow(data)
  stopifnot(d == ncol(data) - 1)
  stopifnot(!is.null(h))

  # Move ahead on kms
  step_ahead <- function(y) {

    # Projected gradient
    if (is_numDeriv) {

      dkde <- numDeriv::grad(func = function(z) {
        DirStats::kde_dir(x = z / sqrt(sum(z^2)), data = data, h = h)
      }, x = y) # TODO: replace with analytical computation
      kde <- DirStats::kde_dir(x = y, data = data, h = h)
      eta <- dkde / kde

    } else {

      eta <- polykde::grad_hess_kde_polysph(x = rbind(y  / sqrt(sum(y^2))),
                                               X = data, d = d, h = h,
                                               norm_grad_hess = TRUE)$grad

    }

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

  # Cut by the number of unique clusters
  if (!is.null(init_clusters) && is_cut){

    labels <- cutree(tree, k = length(unique(init_clusters)))

  }

  modes <- sapply(seq_len(max(labels)), function(i) {
    rotasym::spherical_mean(data = y[labels == i, , drop = FALSE])
  })
  modes <- t(modes)

  labels_rle <- rle(labels)
  antimodes <- c(1, cumsum(labels_rle$lengths))

  # Return end points
  return(list(end_points = y, cluster = labels, modes = modes,
              antimodes = antimodes, paths = paths, tree = tree, h = h,
              eval.points = x, labels_rle = labels_rle))

}

#' @title Compute the bandwidth
#'
#' @description Calculate the bandwidth associated for the given parameters.
#'
#' @param x. A matrix of size \code{c(nx, d + 1)} with the initial points.
#' @param type. The specific way of calculating the bandwidth: hpi linear or rot up.
#' Defaults to \code{"rot_up"}.
#' @export
bw_kms <- function(x, type = c("hpi_linear_s1", "rot_up")[2]) {
  stopifnot(type == "hpi_linear_s1" || type == "rot_up")
  n <- nrow(x)
  d <- ncol(x) - 1
  if (type == "hpi_linear_s1") {

    # Obtain optimal plug-in bandwidth for density derivative estimation,
    # disregarding periodicity
    return(ks::hpi(DirStats::to_rad(x), deriv.order = 1))

  } else {

    # Upscaled EMI bandwidth for derivative
    fit_mix <- DirStats::bic_vmf_mix(data = x, kappa_max = 1e3)
    return(DirStats::bw_dir_emi(data = x, fit_mix = fit_mix)$h_opt *
             n^(1 / (d + 4)) * n^(-1 / (d + 6)))

  }
}


#' @title Plot the density chart of the linear data variable
#'
#' @description Produces a density plot for the given parameters.
#'
#' @param x. A matrix of size \code{c(nx, d + 1)} with the initial points.
#' @param h. Defaults to \code{NULL}. If it is not filled in, the kernel mean shift is
#' calculated with h equal to \code{ks::hpi(x, deriv.order = 1)}.
#' @param init_clusters. Array with the label associated to each observation.
#' Defaults to \code{NULL}.
#' @param tol tolerance for equality of modes. Defaults to \code{1e-1}.
#' @param step. Separation between evaluation points. Defaults to \code{0.01}.
#' @param ylim. Limit to the y coordinate of the plot. Defaults to \code{c(0, 0.35)}.
#' @param is_cut. Boolean value that informs whether the tree can be cut or not.
#' @return Same object that the function \code{kms_dir} returns.
#' Defaults to \code{TRUE}.
#' @export
plot_kde <- function(x, h, tol = 1e-1, init_clusters = NULL, step = 0.01,
                     ylim = c(0, 0.35), is_cut = TRUE) {
  stopifnot(!is.null(x))
  stopifnot(!is.null(h))

  # Convert to radians the scores
  rad <- DirStats::to_rad(x)

  # Hack to have periodicity
  samp_rad <- c(rad - 2 * pi, rad, rad + 2 * pi)

  # Convert to cartesian coordinates (circumference)
  samp_x <- DirStats::to_cir(samp_rad)

  # Generate the evaluation points
  margin <- 1.25
  x_min <- -margin * pi
  x_max <- margin * pi
  eval.points.rad <- seq(x_min, x_max, by = step)

  # Convert them into Cartesian coordinates
  eval.points <- DirStats::to_cir(eval.points.rad)

  # Compute the kernel mean shift
  kms_data <- kms_dir(data = samp_x, x = eval.points, h = h,
                      init_clusters = init_clusters, tol = tol, is_cut = is_cut)

  # Set only the modes between [-margin*pi, margin*pi]
  modes <- DirStats::to_rad(kms_data$modes)
  modes <- c(modes -2 * pi, modes, modes + 2 * pi)
  modes <- modes[modes >= x_min & modes <= x_max]

  unique_modes <- sort(modes)
  unique_modes <- unique_modes[unique_modes >= -pi & unique_modes <= pi]
  total_modes <- length(unique_modes)

  # Calculate the labels for the interval [-pi, pi] (not the space within the margin)
  labels_interval <- kms_data$cluster[eval.points.rad >= -pi & eval.points.rad <= pi]
  labels_int_rle_values <- rle(labels_interval)

  # Getting the antimodes for all the space, counting the margin
  positions_antimodes <- kms_data$antimodes
  labels_rle_values <- kms_data$labels_rle


  # Colors for the plot, different alpha (more or less opaque)
  col <- rainbow(total_modes + 2, alpha = 1)
  col_alpha <- rainbow(total_modes + 2, alpha = 0.15)


  # Plot curve and axes
  plot(x = eval.points.rad, y = DirStats::kde_dir(x = eval.points, data = x, h = h),
       xlim = c(-pi, pi),
       axes = FALSE, xlab = "", ylab = "",
       ylim = ylim, type = "l"
  )
  sdetorus::torusAxis(1)
  axis(2)
  abline(h = 0, lty = 3)


  # Plot modes
  for (i in seq_along(unique_modes)) {

    segments(
      x0 = unique_modes[i], y0 = 0, x1 = unique_modes[i],
      y1 = DirStats::kde_dir(
        data = x, h = h, x = DirStats::to_cir(unique_modes[i])
      ),
      col = col[labels_int_rle_values$values[i]], lwd = 1
    )

  }

  # Plot domains of attraction
  kde <- DirStats::kde_dir(x = eval.points, data = x, h = h)
  for (k in seq_along(labels_rle_values$values)) {

    begin <- positions_antimodes[k]
    end <- positions_antimodes[k + 1]
    polygon(
      x = c(
        eval.points.rad[begin],
        eval.points.rad[begin:end],
        eval.points.rad[end]
      ),
      y = c(0, kde[begin:end], 0),
      col = col_alpha[labels_rle_values$values[k]],
      border = NA
    )

  }

  lines(eval.points.rad, kde)
  abline(v = eval.points.rad[positions_antimodes], lty = 3)

  pi_interval_index <- samp_rad >= - (margin - 0.17) * pi & samp_rad <= (margin - 0.17) * pi
  samp_rad_interv <- samp_rad[pi_interval_index]
  orig_clusters <- rep(init_clusters, 3)[pi_interval_index]

  if (!is.null(init_clusters)) {

    # Plot rug with colors following the initial clustering
    for (i in unique(init_clusters)) {
      rug_col <- viridis::viridis(length(unique(init_clusters)))
      rug(samp_rad_interv[orig_clusters == i], col = rug_col[i], ticksize = 0.03)

    }

  } else {

    # Plot rug without coloring (clustering)
    rug(samp_rad_interv, ticksize = 0.03)

  }
  return(kms_data)
}
