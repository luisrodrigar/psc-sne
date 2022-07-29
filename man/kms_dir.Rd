% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kms.R
\name{kms_dir}
\alias{kms_dir}
\title{Kernel mean shift clustering for directional data}
\usage{
kms_dir(x, data, h = NULL, N = 500, eps = 0.001, tol = 0.1,
  keep_paths = FALSE, show_prog = TRUE)
}
\arguments{
\item{x}{a matrix of size \code{c(nx, d + 1)} with the initial points for
the Euler algorithm.}

\item{data}{a matrix of size \code{c(n, d + 1)} with the sample.}

\item{h}{bandwidth. Chosen automatically if \code{NULL} (default).}

\item{N}{maximum number of iterations. Defaults to \code{500}.}

\item{eps}{convergence tolerance. Defaults to \code{1e-3}.}

\item{tol}{tolerance for equality of modes. Defaults to \code{1e-1}.}

\item{keep_paths}{keep the ascending paths? Defaults to \code{FALSE}.}

\item{show_prog}{display progress?}
}
\value{
A list with the following entries:
\itemize{
  \item \code{end_points}: end points of the Euler algorithm. A matrix of
  the same size as \code{x}.
  \item \code{cluster}: vector giving the cluster labels.
  \item \code{modes}: estimated modes for each cluster (sorted).
  \item \code{paths}: ascension paths, if \code{keep_paths = TRUE}. A list
  of length \code{nx} with matrices of size \code{c(np, d + 1)}, where
  \code{np} is at most \code{N + 1}.
}
}
\description{
Performs kernel mean shift clustering on \eqn{S^d} using
an adapted Euler algorithm and kernel density estimator.
}
\examples{
# Detection of three clusters in S^2
data <- rbind(
  rotasym::r_vMF(n = 50, mu = c(0, 0, 1), kappa = 5),
  rotasym::r_vMF(n = 50, mu = c(0, 0, -1), kappa = 5),
  rotasym::r_vMF(n = 50, mu = c(1, 0, 0), kappa = 5)
)
kms <- kms_dir(x = data, data = data, keep_paths = TRUE)
sd3 <- scatterplot3d::scatterplot3d(data, xlim = c(-1, 1),
                                    ylim = c(-1, 1), zlim = c(-1, 1),
                                    color = kms$cluster + 1, pch = 16,
                                    cex.symbol = 0.5)
for (i in seq_len(nrow(data))) sd3$points3d(kms$paths[[i]], type = "l",
                                            lty = 3)
sd3$points3d(kms$end_points, col = kms$cluster + 1, pch = "*", cex = 2)

# Detection of three clusters in S^1
data <- rbind(
  rotasym::r_vMF(n = 50, mu = c(0, 1), kappa = 5),
  rotasym::r_vMF(n = 50, mu = c(-sqrt(2), -sqrt(2)) / 2, kappa = 5),
  rotasym::r_vMF(n = 50, mu = c(sqrt(2), -sqrt(2)) / 2, kappa = 5)
)
kms <- kms_dir(x = data, data = data, keep_paths = TRUE)
plot(data, col = kms$cluster + 1, pch = 16, xlim = c(-1.5, 1.5),
     ylim = c(-1.5, 1.5))
for (i in seq_len(nrow(data))) {
 l <- seq(0, 1, length.out = nrow(rbind(kms$paths[[i]])))
 lines(sqrt(1 + l) * kms$paths[[i]], lty = 3)
}
points(sqrt(2) * kms$end_points, col = kms$cluster + 1, pch = "*", cex = 2)
kms
}