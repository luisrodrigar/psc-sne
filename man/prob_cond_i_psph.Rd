% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_distribution.R
\name{prob_cond_i_psph}
\alias{prob_cond_i_psph}
\title{Conditional polyspherical Cauchy probability for the
\eqn{i}-th observation}
\usage{
prob_cond_i_psph(x, i, rho_list, d_total_i_psph_cauchy = NULL)
}
\arguments{
\item{x}{an array of size \code{c(n, d + 1, r)} with the polyspherical data,
where \code{n} is the number of observations, \code{d} is the dimension of
each sphere, and \code{r} is the number of spheres.}

\item{i}{corresponds to the \eqn{i}-th observation index.}

\item{rho_list}{rho list of size \code{n} for each \eqn{i}-th observation
that stands for the concentration parameter.}

\item{d_total_i_psph_cauchy}{marginal probability of the \eqn{i}-th
observation. Optional, defaults to \code{NULL}.}
}
\value{
Conditional polyspherical Cauchy probability for all the \eqn{j}-th given the \eqn{i}-th observation.
}
\description{
Calculate the conditional high-dimension probability
for all the \eqn{j}-th observation given the \eqn{i}-th.
}
\examples{
x <- sphunif::r_unif_sph(20, 3, 4)
rho_list_1 <- rep(0.5, 20)
rho_list_2 <- rep(0.9999, 20)
d_total_1_psph <- d_total_psph_cauchy(x, rho_list_1)[1]
d_total_4_psph <- d_total_psph_cauchy(x, rho_list_2)[4]
prob_cond_i_psph(x, 1, rho_list_1, d_total_1_psph)
prob_cond_i_psph(x, 4, rho_list_2, d_total_4_psph)
}
