
#' @title Cluster factory
#'
#' @description Cluster factory for concurrent tasks.
#'
#' @param num_cores number of cores to be available in the cluster.
#' @param outfile text file where print the output, empty means no output.
#' @return Cluster created with the number of cores passing as parameter.
#' @export
#' @examples
#' clusterFactory(2)
#' clusterFactory(2, "log.txt")
#' file.remove("log.txt")
clusterFactory <- function(num_cores, outfile = "") {
  cl <- NULL
  if (tolower(.Platform$OS.type) != "windows") {
    cl <- parallel::makeCluster(spec = num_cores, type = "FORK", outfile = outfile)
  } else {
    cl <- parallel::makeCluster(spec = num_cores, outfile = outfile)
  }
  return(cl)
}

#' @title Reconstruct cosine similarities matrix
#'
#' @description Reconstruct the cosine similarities from vector elements to symmetric matrix.
#'
#' @param cos_sim_vec vector containing all non-repeated values of the cosine similarity matrix (upper triangular matrix).
#' @param n sample size of the original data.
#' @return Cosine similarity matrix formed by the \code{cos_sim_vector}.
#' @export
#' @examples
#' n <- 6
#' x <- sphunif::r_unif_sph(n, 3)
#' cos_sim_vec <- sphunif::Psi_mat(x, scalar_prod = TRUE)
#' reconstruct_cos_sim_mat(cos_sim_vec[, 1], n)
reconstruct_cos_sim_mat <- function(cos_sim_vec, n) {
  cos_sim_mat <- matrix(0, nrow = n, ncol = n)
  cos_sim_mat[upper.tri(cos_sim_mat)] <- cos_sim_vec
  cos_sim_mat <- cos_sim_mat + t(cos_sim_mat)
  diag(cos_sim_mat) <- 1
  return(cos_sim_mat)
}

#' @title Cosine similarity for the poly-sphere
#'
#' @description Calculates the cosine similarity for each sphere of the 3 dimensional array.
#'
#' @inheritParams high_dimension
#' @return An array of size \code{c(n, n, r)} with the cosine similarities of each sphere \eqn{\mathcal{S}^d} from \code{x}.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(100, 3, 3)
#' cosine_polysph(x)
cosine_polysph <- function(x) {
  r <- dim(x)[3]
  n <- nrow(x)
  # applying the cosine function to each third dimension of the array
  cos_sim_vec <- sphunif::Psi_mat(x, scalar_prod = TRUE)
  sapply(1:r, function(k) reconstruct_cos_sim_mat(cos_sim_vec[, k], n), simplify = "array")
}

#' @title Symmetric probabilities
#'
#' @description Calculate the symmetric probabilities of a given conditional polyspherical Cauchy probability matrix.
#'
#' @param P matrix of probabilities \eqn{(P_{i|j})_{ij}} with size \code{c(n, n)} where \code{n} is the number of observations of the original array \code{x}.
#' @return The sum of \code{P} and \code{t(P)} divided by twice the sample size.
#' @export
#' @examples
#' symmetric_probs(matrix(runif(3 * 3), nrow = 3, ncol = 3))
#' symmetric_probs(diag(3))
symmetric_probs <- function(P) {
  n <- nrow(P)
  P <- (P + t(P)) / (2 * n)
  return(P)
}

#' @title Set each object's diagonal of a 3d-array
#'
#' @description Set a specific value to the diagonal of a each matrix of a 3 dimensional array.
#'
#' @inheritParams high_dimension
#' @param k the \code{k}-th sphere from \eqn{(\mathcal{S}^p)^r}.
#' @param val value to set to each array's diagonal of the poly-sphere \code{x}.
#' @return The \code{k}-th matrix of \code{x} with the diagonal set to \code{val}.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(100, 3, 3)
#' diag_3d(x, 1, 0)
#' diag_3d(x, 3, 0)
diag_3d <- function(x, k, val) {
  if (k < 1 || k > dim(x)[3]) {
    stop("The 3rd dimensional index k not valid, must be >= 1 and <= r")
  }
  diag(x[, , k]) <- val
  return(x[, , k])
}

#' @title Radial projection onto the sphere
#'
#' @description Projection of the points onto the sphere of radius 1.
#'
#' @param y matrix with the points in the sphere.
#' @return A matrix with the values of x projected onto the sphere of radius 1.
#' @export
#' @examples
#' y <- rotasym::r_unif_sphere(100, 2)
#' radial_projection(y)
radial_projection <- function(y) {
  # For every observation, apply the formula of the radial projection
  t(sapply(1:nrow(y), function(i) {
    # y_i / |y_i|
    y[i, ] / norm(y[i, ], type = "2")
  }))
}

#' @title Generate optimum evenly separated points
#'
#' @description Generated optimal evenly separated points onto the sphere \eqn{\mathcal{S}^d}.
#'
#' @param n positive integer that defines the size of the sample to generate.
#' @param d size of the low-dimension which defines the sphere \eqn{\mathcal{S}^d}.
#' @return Evenly optimal separated points onto the low-dimension sphere \eqn{\mathcal{S}^d}.
#' @export
#' @examples
#' gen_opt_sphere(100, 1)
#' gen_opt_sphere(250, 2)
gen_opt_sphere <- function(n, d) {
  Y <- NULL
  if (n < 1) {
    stop("n not valid, must be a positive integer")
  }
  # Generate data for the circumference S^1
  if (d == 1) {
    Y <- DirStats::to_cir(seq(0, 2 * pi, l = n + 1)[-(n + 1)]) # 0 == 2 * pi, so we exclude it
  }
  # Generate data for the sphere S^2
  else if (d == 2) {
    Y <- fibonacci_lattice(n)
  } else {
    stop("d not valid, only 1 or 2 are allowed")
  }
  return(Y)
}

#' @title Fibonacci lattice algorithm to generate evenly separated points
#'
#' @description Generated optimal evenly separated points onto the sphere \eqn{\mathcal{S}^2}.
#'
#' @inheritParams gen_opt_sphere
#' @return Evenly optimal separated points onto the low-dimension sphere \eqn{\mathcal{S}^2}.
#' @export
#' @examples
#' fibonacci_lattice(100)
#' fibonacci_lattice(250)
fibonacci_lattice <- function(n) {
  i <- seq(0, n - 1) + 0.5
  phi <- acos(1 - 2 * i / n)
  goldenRatio <- (1 + 5**0.5) / 2
  theta <- 2 * pi * i / goldenRatio
  x <- cos(theta) * sin(phi)
  y <- sin(theta) * sin(phi)
  z <- cos(phi)
  return(cbind(x, y, z))
}
