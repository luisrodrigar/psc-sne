
#' Cluster factory
#'
#' @param num_cores number of cores
#' @param outfile text file where print the output, empty means no output
#' @return the cluster created
#' @examples
#' clusterFactory(4)
#' clusterFactory(7, "log.txt")
clusterFactory <- function(num_cores, outfile = "") {
  cl <- NULL
  if (tolower(.Platform$OS.type) != "windows") {
    cl <- makeCluster(spec = num_cores, type = "FORK", outfile = outfile)
  } else {
    cl <- makeCluster(spec = num_cores, outfile = outfile)
  }
  return(cl)
}

#' Poly-spherical data generation
#'
#' @param n sample size
#' @param p number of dimension in the sphere (S^p)
#' @param r number of spheres (S^p)^r
#' @return data onto the sphere
#' @examples
#' gen_polysphere(100, 1, 3)
#' gen_polysphere(250, 2, 5)
gen_polysphere <- function(n, p, r) {
  polysphere <- array(NA, dim = c(n, p + 1, r))
  for (k in seq_len(r)) {
    polysphere[, , k] <- r_unif_sphere(n, p + 1)
  }
  polysphere
}

#' Generated optimal evenly separated points onto the sphere S^d
#'
#' @param n number of points to be generated
#' @param d reduced dimension
#' @return evenly optimal separated points onto the low-dimension sphere S^d
#' @examples
#' gen_opt_sphere(100, 1)
#' gen_opt_sphere(250, 2)
gen_opt_sphere <- function(n, d) {
  Y <- NULL
  if (n < 1) {
    stop("Parameter n not valid, it has to be positive integer")
  }
  if (d == 1) {
    Y <- DirStats::to_cir(seq(0, 2 * pi, l = n + 1)[-(n + 1)]) # 0 == 2 * pi, so we exclude it
  } else if (d == 2) {
    Y <- X_Fib(n)
  } else {
    stop("Parameter d not valid, only 1 or 2")
  }
  return(Y)
}

X_Fib <- function(N) {
  phi_inv <- (1 - sqrt(5)) / 2
  N <- N / 2
  i <- (-N:N)[-N]
  lat <- asin(2 * i / (2 * N + 1))
  lon <- (2 * pi * i * phi_inv) %% (2 * pi)
  X <- cbind(cbind(cos(lon), sin(lon)) * cos(lat), sin(lat))
  return(X)
}
