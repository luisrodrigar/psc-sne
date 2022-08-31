

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
#'   \item{lat}{latitude coordinate in decimal format}
#'   \item{lon}{longitude coordinate in decimal format}
#'   \item{time}{time in the following format yyyy-MM-dd hh:mm:ss}
#'   \item{theta}{angle in radians of the sea current vector}
#'   \item{speed}{module of the sea current vector}
#' }
#' @details
#' The data object is created thanks to the code in
#' \url{https://github.com/luisrodrigar/psc-sne/blob/main/data-raw/strait-juan-fuca/data-acquisition-juan-fuca.R}
#' @source \url{https://github.com/luisrodrigar/psc-sne/blob/main/data/jdf.rda}
#' @references
#' TODO
#' @examples
#' # Load data
#' data("jdf")
#'
#' # TODO
"jdf"


#' @title Small RNA dataset
#'
#' @description TODO
#'
#' @docType data
#' @format A data frame with TODO rows and TODO variables:
#' \describe{
#'   \item{TODO}{TODO}
#'   \item{TODO}{TODO}
#' }
#' @details
#' TODO
#' @source TODO
#' @references
#' TODO
#' @examples
#' # Load data
#' data("smallrna")
#'
#' # TODO
"smallrna"
