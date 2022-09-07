

#' @title Juan de Fuca currents
#'
#' @description Sea currents defined by their angle and speed produced
#' in some geographic coordinates in Juan de Fuca strait area from
#' 2020-06-01 to 2022-07-01, last day excluded.
#'
#' @docType data
#' @format A data frame with 638,400 observations (rows) and 5 variables
#' (columns):
#' \describe{
#'   \item{lat}{latitude coordinate in decimal format.}
#'   \item{lon}{longitude coordinate in decimal format.}
#'   \item{time}{time in the following format yyyy-MM-dd hh:mm:ss.}
#'   \item{theta}{angle in radians of the sea current vector.}
#'   \item{speed}{module of the sea current vector.}
#' }
#' @details
#' The data object is created with the code in
#' \url{https://github.com/luisrodrigar/psc-sne/blob/main/data-raw/strait-juan-fuca/data-acquisition-juan-fuca.R}
#' @source \url{https://github.com/luisrodrigar/psc-sne/blob/main/data/jdf.rda}
#' @references
#' TODO
#' @examples
#' # Load data
#' data("jdf")
#'
#' # Plot vector field
#' # TODO
"jdf"


#' @title Small RNA dataset
#'
#' @description "Small RNA dataset". Among others, used in Section 5.2 in
#' Zoubouloglou et al. (2022) and references therein.
#'
#' @docType data
#' @format A data frame with 190 rows and 10 variables:
#' \describe{
#'   \item{angles}{matrix of 7 dihedral angles.}
#'   \item{clusters}{vector of cluster labels.}
#'   \item{torsion}{matrix of 2 torsion angles.}
#' }
#' @details
#' Clusters are defined from the torsion angles, see Section 5.2 in
#' Zoubouloglou et al. (2022) and references therein.
#' @source Dataset put together by Duarte and Pyle (1998) and updated by Wadley
#' et al. (2007). Previously used by Eltzner et al. (2018).
#' @references
#' Duarte, C. M. and Pyle, A. M. (1998). Stepping through an RNA structure: A
#' novel approach to conformational analysis. \emph{Journal of Molecular
#' Biology}, 284(5):1465--1478. \doi{10.1006/jmbi.1998.2233}.
#'
#' Eltzner, B., Huckemann, S., and Mardia, K. V. (2018). Torus principal
#' component analysis with applications to RNA structure. \emph{The Annals of
#' Applied Statistics}, 12(2):1332--1359. \doi{10.1214/17-AOAS1115}.
#'
#' Wadley, L. M., Keating, K. S., Duarte, C. M., and Pyle, A. M. (2007).
#' Evaluating and learning from RNA pseudotorsional space: Quantitative
#' validation of a reduced representation for RNA structure. \emph{Journal of
#' Molecular Biology}, 372(4):942--957. \doi{10.1016/j.jmb.2007.06.058}.
#'
#' Zoubouloglou, P., García-Portugués, E., and Marron, J. S. (2022). Scaled
#' torus principal component analysis. \emph{Journal of Computational and
#' Graphical Statistics}. \doi{10.1080/10618600.2022.2119985}.
#' @examples
#' # Load data
#' data("smallrna")
#'
#' # Clusters
#' pairs(smallrna$angles, col = smallrna$clusters, pch = 16)
"smallrna"
