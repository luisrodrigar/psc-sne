% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{index_upper_trian}
\alias{index_upper_trian}
\title{Vector index of an upper triangular matrix without diagonal}
\usage{
index_upper_trian(i, j, n)
}
\arguments{
\item{i}{corresponds to the \eqn{i}-th observation index.}

\item{j}{corresponds to the \eqn{j}-th observation index.}

\item{n}{the sample size of the symmetric matrix of size \code{c(n, n)}.}
}
\description{
Obtain the index of an element within a vector
for the equivalent in the upper triangular matrix without diagonal
}
\examples{
mat <- rbind(c(0, 1, 2), c(0, 0, 1), c(0, 0, 0))
vec <- c(1, 2, 1)
index <- index_upper_trian(1, 2, 3)
all.equal(vec[index], mat[1, 2])
}
\keyword{internal}
