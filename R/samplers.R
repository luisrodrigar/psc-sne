

#' @title Samplers of one-dimensional modes of variation for polyspherical data
#'
#' @description Functions for sampling data on \eqn{(S^d)^r} for \eqn{d=1,2}
#' using one-dimensional modes of variation.
#'
#' @param n sample size.
#' @param r number of spheres in the polysphere \eqn{(S^d)^r}.
#' @param alpha a vector of size \code{r} valued in \eqn{[-\pi,\pi)} with the
#' initial angles for the linear trend. Chosen at random by default.
#' @param k a vector of size \code{r} with the \bold{integer} slopes defining
#' the angular linear trend. Chosen at random by default.
#' @param sigma standard deviation of the noise about the one-dimensional mode
#' of variation. Defaults to \code{0.25}.
#' @param angles return angles in \eqn{[-\pi, \pi)}? Defaults to \code{FALSE}.
#' @param c \href{https://en.wikipedia.org/wiki/Clélie}{Clélie curve}
#' parameter, changing the spiral wrappings. Defaults to \code{1}.
#' @param t latitude, with respect to \code{Theta}, of the small circle.
#' Defaults to \code{0} (equator).
#' @param Theta a matrix of size \code{c(3, r)} giving the north poles for
#' \eqn{S^2}. Useful for rotating the sample. Chosen at random by default.
#' @param spiral consider a spiral (or, more precisely, a
#' \href{https://en.wikipedia.org/wiki/Clélie}{Clélie curve}) instead of
#' a small circle? Defaults to \code{FALSE}.
#' @return
#' An array of size \code{c(n, d, r)} with samples on \eqn{(S^d)^r}. If
#' \code{angles = TRUE} for \code{r_path_s1r}, then a matrix of size
#' \code{c(n ,r)} with angles is returned.
#' @examples
#' # Straight trends on (S^1)^2
#' n <- 200
#' samp_1 <- r_path_s1r(n = n, r = 2, k = c(1, 2), angles = TRUE)
#' plot(samp_1, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n),
#'      axes = FALSE, xlab = "", ylab = "", pch = 16)
#' sdetorus::torusAxis()
#'
#' # Straight trends on (S^1)^3
#' n <- 200
#' samp_2 <- r_path_s1r(n = n, r = 3, angles = TRUE)
#' pairs(samp_2, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n),
#'       pch = 16)
#' sdetorus::torusAxis()
#'
#' # Small-circle trends on (S^2)^2
#' n <- 200
#' samp_3 <- r_path_s2r(n = n, r = 2, sigma = 0.1)
#' old_par <- par(mfrow = c(1, 2))
#' scatterplot3d::scatterplot3d(
#'   samp_3[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' scatterplot3d::scatterplot3d(
#'   samp_3[, , 2], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' par(old_par)
#'
#' # Spiral trends on (S^2)^2
#' n <- 200
#' samp_4 <- r_path_s2r(n = n, r = 2, c = 3, spiral = TRUE, sigma = 0.01)
#' old_par <- par(mfrow = c(1, 2))
#' scatterplot3d::scatterplot3d(
#'   samp_4[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' scatterplot3d::scatterplot3d(
#'   samp_4[, , 2], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' par(old_par)
#' @name samplers_one


#' @rdname samplers_one
#' @export
r_path_s1r <- function(n, r, alpha = runif(r, -pi, pi),
                       k = sample(-2:2, size = r, replace = TRUE),
                       sigma = 0.25, angles = FALSE) {

  # Angular trend
  t <- sort(runif(n))
  theta <- t(alpha + (2 * pi * k) %o% t)

  # Add noise
  eps <- matrix(rnorm(r * n, sd = sigma), nrow = n, ncol = r)
  theta <- sdetorus::toPiInt(theta + eps)

  # Cartesianize
  if (angles) {

    return(theta)

  } else {

    return(sphunif::Theta_to_X(theta))

  }

}


#' @rdname samplers_one
#' @export
r_path_s2r <- function(n, r, t = 0, c = 1,
                       Theta = t(rotasym::r_unif_sphere(n = r, p = 3)),
                       sigma = 0.25, spiral = FALSE) {

  samp <- array(dim = c(n, 3, r))
  if (spiral) {

    # Sample common latitude
    ph <- sort(runif(n, min = 0, max = pi))

    # Loop on the S^2's
    for (j in 1:r) {

      # Unrotated sample
      eps <- rnorm(n, sd = sigma)
      samp[, , j] <- DirStats::to_sph(th = c * ph + eps, ph = ph)

      # Rotated sample
      samp[, , j] <- samp[, , j] %*%
        cbind(Theta[, j], rotasym::Gamma_theta(theta = Theta[, j]))

    }

  } else {

    # Sample common and unrotated longitudes
    th <- sort(runif(n, min = -pi, max = pi))
    U <- DirStats::to_cir(th = th)

    # Loop on the S^2's
    for (j in 1:r) {

      # Use the moment estimators of a Beta(shape1, shape2) to get the
      # parameters yielding
      # m = E[Beta(shape1, shape2)] = (t + 1) / 2 and
      # v = Var[Beta(shape1, shape2)] = sigma^2 / 4
      # (https://en.wikipedia.org/wiki/Beta_distribution#Two_unknown_parameters)
      # In this way,
      # E[2 * Beta(shape1, shape2) - 1] = t and
      # Var[2 * Beta(shape1, shape2) - 1] = sigma^2.
      m <- (t + 1) / 2
      v <- 0.25 * sigma^2
      stopifnot(v < m * (1 - m))
      shape1 <- m * (m * (1 - m) / v - 1)
      shape2 <- (1 - m) * (m * (1 - m) / v - 1)

      # Sample latitudes centered at t and with standard deviation sigma
      V <- 2 * rbeta(n = n, shape1 = shape1, shape2 = shape2) - 1

      # Tangent-normal decompositon about Theta[, j]
      samp[, , j] <- V * matrix(Theta[, j], nrow = n, ncol = 3, byrow = TRUE) +
        sqrt(1 - V^2) * U %*% t(rotasym::Gamma_theta(theta = Theta[, j]))

    }

  }

  return(samp)

}


#' @title Sample correlation matrices with block structure
#'
#' @description Sample a zero-mean multivariate normal
#' \eqn{N_{gp}(0, \Sigma)} with a diagonal block matrix \eqn{\Sigma}
#' partitioned into \eqn{g} blocks with \eqn{p} variables. Each block is
#' constructed as a \code{\link[=stats]{toeplitz}} correlation matrix.
#'
#' @inheritParams samplers_one
#' @param g number of groups of variables that are uncorrelated between them.
#' @param p number of variables on each group.
#' @param rho a vector of size \code{g} with the correlations determining the
#' \code{\link[=stats]{toeplitz}} correlation matrices of each group.
#' @return A matrix of size \code{c(n, g * p)} with the sample.
#' @examples
#' # Visualize some features
#' n <- 200
#' g <- 10
#' p <- 10
#' x <- r_block(n = n, g = g, p = p)
#' pairs(x[, c(1:2, p + 1:2)])
#'
#' # Standardize variables -- now the vectors of observations for each variable
#' # (the columns) live on \sqrt{n - 1} * S^{n - 1}!
#' x_sca <- scale(x)
#' colSums(x_sca^2)
#'
#' # Make the features live on S^{n - 1}
#' x_sca <- scale(x) / sqrt(n - 1)
#'
#' # Transpose matrix (features become observations)
#' feat <- t(x_sca)
#'
#' # Run psc_sne() with colors being the groups of variables
#' # dim(feat) <- c(dim(feat), 1)
#' # cols <- rep(1:g, each = p)
#' # psc_sne(X = feat, d = 1, colors = cols)
#' @export
r_block <- function(n, g = 5, p = 20, rho = rep(c(-0.9, 0.9), times = g)[1:g]) {

  # Check rho
  if (length(rho) != g) {

    stop("rho must have length g")

  }

  # Blockwise covariance matrix
  big_p <- g * p
  Sigma <- matrix(0, nrow = big_p, ncol = big_p)
  for (i in 1:g) {

    ind_i <- (1 + p * (i - 1)):(p * i)
    Sigma[ind_i, ind_i] <- toeplitz(rho[i]^(0:(p - 1)))

  }

  # Sample N_{g * p}(0, Sigma)
  return(mvtnorm::rmvnorm(n = n, mean = rep(0, big_p), sigma = Sigma))

}


#' #' @title Sample uniform data on the polysphere
#' #' @export
#' r_unif_psph <- function(n, d, r, angles = FALSE) {
#'
#'   if (angles && d == 1) {
#'
#'     samp <- sdetorus::toPiInt(sphunif::r_unif_cir(n = n))
#'
#'   } else {
#'
#'     samp <- sphunif::r_unif_sph(n = n, p = d + 1, M = r)
#'
#'   }
#'   return(samp)
#'
#' }
