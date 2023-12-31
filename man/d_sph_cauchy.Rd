% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_distribution.R
\name{d_sph_cauchy}
\alias{d_sph_cauchy}
\title{High-dimension polyspherical Cauchy density}
\usage{
d_sph_cauchy(x, i, j, rho, k, p)
}
\arguments{
\item{x}{an array of size \code{c(n, d + 1, r)} with the polyspherical data,
where \code{n} is the number of observations, \code{d} is the dimension of
each sphere, and \code{r} is the number of spheres.}

\item{i}{corresponds to the \eqn{i}-th observation index.}

\item{j}{corresponds to the \eqn{j}-th observation index.}

\item{rho}{concentration parameter must be in [0, 1).}

\item{k}{corresponds to the \eqn{k}-th sphere index.}

\item{p}{dimension of the sphere, \eqn{\mathcal{S}^p}.}
}
\value{
Spherical Cauchy density value given the parameters.
}
\description{
Calculate the high-dimension spherical Cauchy density function.
}
