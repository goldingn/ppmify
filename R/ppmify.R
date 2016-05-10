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
#' @import sp

NULL


#' @name ppmify
#' @title create a ppm object from point data
#' @export
#' @param coords a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in the PPM analysis. If a matrix or
#'   dataframe, the first column should give the horizontal
#'   (x/longitude/easting) coordinates, and the second column the vertical
#'   (y/latitude/northing) coordinates.
#' @param area an optional extent, SpatialPolygons* or Raster* object giving the
#'   area over which to model the point process. If ignored, a rectangle
#'   defining the extent of \code{coords} will be used instead.
#' @param covariates an optional Raster* object containing covariates for
#'   modelling the point process
#' @param method the method for selecting quadrature points. This will either
#'   generate a set of integration points with appropriate weights, or count the
#'   number of points falling in each cell (if \code{method = 'count'}). See
#'   \code{Details} for information on the available approaches.
#' @param density the number of integration points required per square kilometre
#'   (ignored if \code{method = 'count'})
#' @description create a \code{ppm} object containing the information needed to
#'   fit a Poisson point process model using Poisson regression modelling
#'   software.
#' @details <integration details to be added>
#' @return an object of classes \code{ppm} and \code{data.frame}
#'
#' @examples
#' # generate some fake point data
#' r <- raster(system.file("external/test.grd", package="raster"))
#' pts <- sampleRandom(r, 100, xy = TRUE)[, 1:2]
#' plot(r, col = grey(0.8))
#' points(pts, pch = 16, cex = 0.5)
#'
#' # generate ppm data
#' ppm <- ppmify(pts, area = r, covariates = r)
#'
#' # fit a model
#' m <- glm(points ~ offset(offset) + test,
#'          data = ppm,
#'          family = poisson)
#'
#' # predict to a raster, remembering to set the offset value
#' p <- predict(r, m,
#'  const = data.frame(offset = 0),
#'  type = 'response')
#'
#' # plot results (prediction is points per square km)
#' plot(p)
#' points(ppm[ppm$points == 1, c('x', 'y')], pch = 16, cex = 0.5)


ppmify <- function (coords,
                    area = NULL,
                    covariates = NULL,
                    method = c('grid', 'count'),
                    density = 10) {

  # get the requested method
  method <- match.arg(method)

  # convert coordinates into a dataframe
  coords <- coords2df(coords)

  # if an area isn't provided, generate a bounding box
  if (is.null(area)) area <- defaultArea(coords)

  # generate integration points
  int <- switch(method,
                grid = grid(area, density),
                count = NULL)

  npts <- nrow(coords)
  nint <- nrow(int)

  # set up dataframe
  ppm <- data.frame(points = rep(1:0, c(npts, nint)),
                    x = c(coords$x, int$x),
                    y = c(coords$y, int$y),
                    offset = c(rep(0, npts),
                               log(int$weight)))

  # add covariates if they are provided
  if (!is.null(covariates)) {
    # extract covariates
    covs <- data.frame(extract(covariates, ppm[, c('x', 'y')]))
    names(covs) <- names(covariates)

    # make sure there aren't naming conflicts
    names(covs) <- ifelse(names(covs) %in% names(ppm),
                          paste0(names(covs), '.1'),
                          names(covs))

    # add to ppm
    ppm <- cbind(ppm, covs)

  }

  # define the class
  class(ppm) <- c(class(ppm), 'ppm')

  return (ppm)

}

# convert coordinates (in whatever format they arrive in) into a dataframe
# with the expected columns
coords2df <- function (coords) {

  # check object classes
  expectClasses(coords,
                c('matrix',
                  'data.frame',
                  'SpatialPoints',
                  'SpatialPointsDataFrame'),
                name = 'coords')

  if (is.matrix(coords) | is.data.frame(coords)) {

    # for matrices/dataframes, make sure there are only two columns
    if (ncol(coords) != 2) {
      stop (sprintf('coords should have only two columns, giving the horizontal (x/longitude) then vertical (y/latitude) coordinates. The object passed had %i columns',
                    NCOL(coords)))
    }

    # otherwise, coerce into a data.frame and rename the columns
    df <- data.frame(coords)

  } else {
    # otherwise, for SpatialPoints* objects, just grab the coordinates
    df <- data.frame(coords@coords)
  }

  # set column names
  colnames(df) <- c('x', 'y')

  return (df)

}

# grid method for generating integration points. area is the object passed to
# ppmify and density is the required number of integration points per square
# kilometre
grid <- function (area, density) {

  # generate the grid as a raster, find the centroids and get integration
  # weights by a zonal operation

  # get extent of area as a vector
  ext <- as.vector(extent(area))

  # generate an evenly sized raster grid
  # if the area is in lat/long this is trickier
  if (isLonLat(area)) {

    stop ('lat/long areas not yet supported')

    # get cell areas and total extent area

    # work out implied height from ratio

    # put down the correct number of cells

    # extract their area

  } else {

    # otherwise, it should be in metres and evenly sized cells, so easier

    # get dimensions in km
    dim_area <- abs(diff(ext)[c(1, 3)]) * 0.1 ^ 3

    # get number of grid cells in each direction
    ncells <- round(dim_area * density)

    # build a raster grid matching this
    grid <- raster(extent(area),
                   nrows = ncells[1],
                   ncols = ncells[2])

    grid <- setValues(grid, 1:ncell(grid))

  }

  # extract the surface area of area falling in each cell
  grid <- gridArea(grid, area)

  # get grid coordinates and integration weights
  coords <- data.frame(xyFromCell(grid, 1:ncell(grid)))
  coords$weight <- getValues(grid)

  # remove those with 0 weight
  coords <- coords[coords$weight > 0, ]

  return (coords)

}

# create a default area (bounding box of coordinates)
defaultArea <- function (coords) {
  # get limits
  xlim <- range(coords$x)
  ylim <- range(coords$y)

  # add on 10%
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.1
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.1

  # make a SpatialPolygons object
  p <- Polygon(cbind(x = xlim[c(1, 1, 2, 2)],
                     y = ylim[c(1, 2, 2, 1)]))
  ps <- Polygons(list(p), 1)
  sp <- SpatialPolygons(list(ps))

  return (sp)

}
