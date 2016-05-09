#' @title ppmify package
#' @name ppmify-package
#' @description ppmify is a micropackage to help fitting Poisson point process
#'   models to point data, using one of the very many types Poisson regression
#'   models available in R. The main function is \code{ppmify}, which takes some
#'   coordinates, an area of interest and covariates and produces a \code{ppm}
#'   object, which can be used in place of a dataframe in a Poisson regression
#'   model. See \code{\link{ppmify}} for details.
#' @docType package
#' @import raster

NULL


#' @name ppmify
#' @title create a ppm object from point data
#' @param coords a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in the PPM analysis
#' @param area an optional extent, SpatialPolygons* or Raster* object giving the
#'   area over which to model the point process. If ignored, a rectangle
#'   defining the extent of \code{coords} will be used instead.
#' @param covariates an optional Raster* object containing covariates for
#'   modelling the poit process
#' @param method the method for defining integration area. This will either
#'   generate a set of integration points with appropriate weights, or count the
#'   number of points falling in each cell (if \code{method = 'full'}). See \code{Details} for information on the available approaches.
#'
#' @details <integration details to be added>
#'
#'

ppmify <- function (coords,
                    area = NULL,
                    covariates = NULL,
                    method = c('grid', 'count')) {
  # some function
  ppm <- data.frame(y = NA)


  # define the classes
  class(ppm) <- c(class(ppm), 'ppm')
}
