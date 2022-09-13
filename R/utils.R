
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
#' @keywords internal
clusterFactory <- function(num_cores, outfile = "") {
  cl <- NULL
  if (tolower(.Platform$OS.type) != "windows") {
    cl <- parallel::makeCluster(spec = num_cores, type = "FORK",
                                outfile = outfile)
  } else {
    cl <- parallel::makeCluster(spec = num_cores, outfile = outfile)
  }
  return(cl)
}

#' @title Reconstruct the symmetric matrix without diagonal
#'
#' @description Convert the vector into a symmetric matrix with diagonal
#' equal to a specific value.
#'
#' @param vec vector containing all non-repeated values of the cosine
#' similarity matrix (upper triangular matrix).
#' @param n sample size of the original data.
#' @param diag_value value for the diagonal of the resultant matrix.
#' @return Matrix formed by the \code{vec}.
#' @export
#' @examples
#' n <- 6
#' x <- sphunif::r_unif_sph(n, 3)
#' cos_sim_vec <- drop(sphunif::Psi_mat(x, scalar_prod = TRUE))
#' vec2matrix(cos_sim_vec, n, diag_value = 1)
#' @keywords internal
vec2matrix <- function(vec, n, diag_value) {
  # Create the matrix of dimension n x n
  mat <- matrix(0, nrow = n, ncol = n)
  # Define the elements of the upper triangular matrix as the vector
  mat[upper.tri(mat)] <- vec
  # Convert from a upper triangular matrix into a symmetric one
  mat <- mat + t(mat)
  # Set the elements of the diagonal
  diag(mat) <- diag_value
  return(mat)
}

#' @title Vector index of an upper triangular matrix withou diagonal
#'
#' @description Obtain the index of an element within a vector
#' for the equivalent in the upper triangular matrix without diagonal
#'
#' @inheritParams d_sph_cauchy
#' @param n the sample size of the symmetric matrix of size \code{c(n, n)}.
#' @export
#' @examples
#' mat <- rbind(c(0, 1, 2), c(0, 0, 1), c(0, 0, 0))
#' vec <- c(1, 2, 1)
#' index <- index_upper_trian(1, 2, 3)
#' all.equal(vec[index], mat[1, 2])
#' @keywords internal
index_upper_trian <- function(i, j, n) {
  if (i < 1 || i >= n) {

    stop("i must be within [1, n-1]")

  }

  if (j <= i || j > n) {

    stop("j must be within [i+1, n]")

  }

  return(i + floor(((j - 2) * (j - 1)) / 2))
}

#' @title Cosine similarities vector for the i-th observation
#'
#' @description Obtain the index of an element within a vector
#' for the equivalent in the upper triangular matrix without diagonal
#'
#' @inheritParams high_dimension
#' @inheritParams d_sph_cauchy
#' @param r the number of spheres in the polysphere
#' @param n the sample size.
#' @export
#' @examples
#' n <- 6
#' r <- 6
#' x <- sphunif::r_unif_sph(n, 3, r)
#' cos_sim_psh <- sphunif::Psi_mat(x, scalar_prod = TRUE)
#' cos_sim_i(cos_sim_psh, 4, r, n)
#' @keywords internal
cos_sim_i <- function(cos_sim_psh, i, r, n) {
  sapply(seq_len(r), function(k) {
    sapply(seq_len(n), function(j) {
      if (j == i) {

        # Cosine similarity when the i-th and j-th observation are the same
        return(1)

      }

      if (i > j) {

        # Since it is an upper triangular matrix what represents the vector
        # If the i-th observation is greater than the j-th observation
        # this means that the element is located in the lower triangular matrix
        # Then we swap both indexes
        return(cos_sim_psh[index_upper_trian(j, i, n), k])

      }

      # Get the cosine similarity for the i-th with respect the j-th obs
      return(cos_sim_psh[index_upper_trian(i, j, n), k])
    })
  })
}


#' @title Cosine similarity for the polysphere
#'
#' @description Calculates the cosine similarity for each sphere of
#' the 3 dimensional array.
#'
#' @inheritParams high_dimension
#' @return An array of size \code{c(n, n, r)} with the cosine similarities of
#' each sphere \eqn{\mathcal{S}^d} from \code{x}.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(100, 3, 3)
#' cosine_polysph(x)
#' @keywords internal
cosine_polysph <- function(x) {
  r <- dim(x)[3]
  n <- nrow(x)
  # applying the cosine function to each third dimension of the array
  cos_sim_vec <- sphunif::Psi_mat(x, scalar_prod = TRUE)
  sapply(1:r, function(k) vec2matrix(drop(cos_sim_vec[, k]), n, 1),
         simplify = "array")
}

#' @title Symmetric probabilities
#'
#' @description Calculate the symmetric probabilities of a given conditional
#' polyspherical Cauchy probability matrix.
#'
#' @param P matrix of probabilities \eqn{(p_{i|j})_{ij}} with size
#' \code{c(n, n)} where \code{n} is the number of observations of the original
#' array \code{x}.
#' @return The sum of \code{P} and \code{t(P)} divided by twice the sample size.
#' @export
#' @examples
#' symmetric_probs(matrix(runif(3 * 3), nrow = 3, ncol = 3))
#' symmetric_probs(diag(3))
#' @keywords internal
symmetric_probs <- function(P) {
  n <- nrow(P)
  P <- (P + t(P)) / (2 * n)
  return(P)
}

#' @title Set each object's diagonal of a 3d-array
#'
#' @description Set a specific value to the diagonal of a each matrix of a
#' 3 dimensional array.
#'
#' @inheritParams high_dimension
#' @param k the \code{k}-th sphere from \eqn{(\mathcal{S}^p)^r}.
#' @param val value to set to each array's diagonal of the polysphere \code{x}.
#' @return The \code{k}-th matrix of \code{x} with the diagonal set to
#' \code{val}.
#' @export
#' @examples
#' x <- sphunif::r_unif_sph(100, 3, 3)
#' diag_3d(x, 1, 0)
#' diag_3d(x, 3, 0)
#' @keywords internal
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

  return(y / sqrt(rowSums(rbind(y^2))))

}

#' @title Generate optimum evenly separated points
#'
#' @description Generated optimal evenly separated points on
#' \eqn{\mathcal{S}^1} or \eqn{\mathcal{S}^2}.
#'
#' @param n positive integer with the size of the grid to generate.
#' @param d size of the low-dimension which defines the sphere
#' \eqn{\mathcal{S}^d}, must be 1 or 2.
#' @return For \eqn{\mathcal{S}^1}, evenly optimal separated points.
#' For \eqn{\mathcal{S}^2}, Fibonacci lattice is applied to generate points.
#' \eqn{\mathcal{S}^d}.
#' @export
#' @examples
#' grid_sphere(100, 1)
#' grid_sphere(250, 2)
grid_sphere <- function(n, d) {
  Y <- NULL
  if (n < 1) {
    stop("n not valid, must be a positive integer")
  }
  if (d == 1) {
    # Generate data for the circumference S^1
    # 0 == 2 * pi, so we exclude it
    Y <- DirStats::to_cir(seq(0, 2 * pi, l = n + 1)[-(n + 1)])
  } else if (d == 2) {
    # Generate data for the sphere S^2
    i <- seq(0, n - 1) + 0.5
    phi <- acos(1 - 2 * i / n)
    goldenRatio <- (1 + 5**0.5) / 2
    theta <- 2 * pi * i / goldenRatio
    x <- cos(theta) * sin(phi)
    y <- sin(theta) * sin(phi)
    z <- cos(phi)
    Y <- cbind(x, y, z)
  } else {
    stop("d not valid, only 1 or 2 are allowed")
  }
  return(Y)
}
