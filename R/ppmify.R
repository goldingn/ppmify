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
#' m <- glm(points ~ test + offset(log(weights)),
#'          data = ppm,
#'          family = poisson)
#'
#' # predict to a raster, remembering to set the offset value
#' p <- predict(r, m, type = 'response', const = data.frame(weights = 1))
#'
#' # plot results (prediction is in points per square km)
#' plot(p)
#' points(ppm[ppm$points == 1, c('x', 'y')], pch = 16, cex = 0.5)


ppmify <- function (coords,
                    area = NULL,
                    years = NULL,
                    covariates = NULL,
                    bias_points = NULL,
                    method = c('grid', 'count', 'bias', 'biasRaster'),
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
                count = NULL,
                bias = bias(area, bias_points),
                biasRaster = bias(area, bias_points, offset.raster=TRUE))

  npts <- nrow(coords)
  nint <- nrow(int)

  eps <- sqrt(.Machine$double.eps)

  # set up dataframe
  if (is.null(years)){
    
    ppm <- data.frame(points = rep(1:0, c(npts, nint)),
                      x = c(coords$x, int$x),
                      y = c(coords$y, int$y),
                      weights = c(rep(eps, npts),
                                  int$weight))  
    
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
    
  } else {
    
    ppm <- data.frame(points = rep(1:0, c(npts, nint)),
                      x = c(coords$x, int$x),
                      y = c(coords$y, int$y),
                      year = c(years, sample(years, nint, replace=TRUE)),
                      weights = c(rep(eps, npts),
                                  int$weight))
    
    # add covariates if they are provided
    if (!is.null(covariates)) {
      
      # get index of non temporal covariates
      nt_idx <- which(!(substr(names(covariates),
                               nchar(names(covariates)) - 3,
                               nchar(names(covariates))) %in% c(1900:2100)))
      
      if (length(nt_idx)==0) {
        
        stop('last 4 characters of yearly covariate names must be the year')
      }
      
      # extract non temporal covariates
      nt_covariates <- covariates[[nt_idx]]
      
      # extract temporal covariates
      t_covariates <- covariates[[-nt_idx]]
      
      # find max and min years for temporal covariates
      min_year <- as.numeric(min(substr(names(t_covariates),
                                        nchar(names(t_covariates)) - 3,
                                        nchar(names(t_covariates)))))
      
      max_year <- as.numeric(max(substr(names(t_covariates),
                                        nchar(names(t_covariates)) - 3,
                                        nchar(names(t_covariates)))))
      
      # make list of general names for temporal covs
      general_names <- unique(substr(names(t_covariates), 1, nchar(names(t_covariates))-5))
      
      # get index of temporal covs of max year
      max_idx <- which (substr(names(t_covariates),
                               nchar(names(t_covariates)) - 3,
                               nchar(names(t_covariates)))==max_year)
      
      # create stack of most contemporary covariates
      contemp_covariates <- addLayer(nt_covariates, t_covariates[[max_idx]])
      
      # create empty dataframe with columns for covariates and rows for points
      # extract covariates
      covs <- data.frame(extract(contemp_covariates, ppm[, c('x', 'y')]))
      names(covs) <- c(names(nt_covariates), general_names)
      
      # make empty dataframe to populate later
      covs[] <- NA
      
      # loop through extracting covs for all years
      for (year in min_year:max_year) {
        
        # get years for all presence and pseudo data
        dat_years <- ppm$year
        
        # truncate years to the ones we have covariates for
        dat_years <- pmax(dat_years, min_year)
        dat_years <- pmin(dat_years, max_year)
        
        # index for this year's data
        idx <- which(dat_years == year)
        
        # get index of covs for this year
        idx2 <- which (substr(names(t_covariates),
                              nchar(names(t_covariates)) - 3,
                              nchar(names(t_covariates)))==year)
        
        # get covs for this year
        covs_year <- t_covariates[[idx2]]
        
        # add general names
        names(covs_year) <- general_names
        
        # add nontemporal covariates
        covs_year <- addLayer(nt_covariates, covs_year)
        
        # extract data
        covs_year_extract <- extract(covs_year, ppm[idx, c('x', 'y')])
        
        # check they're all there
        stopifnot(all(colnames(covs) %in% colnames(covs_year_extract)))
        
        # match up the column names so they're in the right order
        match <- match(colnames(covs), colnames(covs_year_extract))
        
        # extract covariates for all points
        covs[idx, ] <- covs_year_extract[, match]
        
        
      } # close years loop 
      
      # add to ppm
      ppm <- cbind(ppm, covs)
      
    } # close !is.null(covariates) statement
  } # close !is.null(years) statement  

  # define the class
  class(ppm) <- c(class(ppm), 'ppm')

  # clean out missing data & warn user
  ppm <- ppmClean(ppm)

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

  # get dimensions in km
  dim_area <- extentDim(ext, lonlat = isLonLat(area))

  # get number of grid cells in each direction
  ncells <- round(dim_area * density)

  # build a raster grid matching this
  grid <- raster(extent(area),
                 nrows = ncells[2],
                 ncols = ncells[1],
                 crs = crs(area))

  # Give cells their number as their value
  grid <- setValues(grid, 1:ncell(grid))

  # extract the surface area of area falling in each cell
  grid <- gridArea(grid, area)

  # get grid coordinates and integration weights
  coords <- data.frame(xyFromCell(grid, 1:ncell(grid)))
  coords$weight <- getValues(grid)

  # remove those with 0 weight
  coords <- coords[coords$weight > 0, ]

  return (coords)

}

bias <- function (area, bias_points, offset.raster=FALSE){
  
  # calculate weights for integration points when using pseudo-observations as integration points
  # weight = total study area / number of integration/quadrature points
  
  if (offset.raster){ # if offset raster is also provided, calculate sum of raster values
    
    a <- cellStats(area, sum, na.rm=TRUE)  
    
  } else {# calculate total study area
    
    a <- cellStats(area(area), sum, na.rm=TRUE)
    
  }
  # calculate number of integration points
  n <- nrow(bias_points)
  
  # calculate weights
  weight <- a/n
  
  # convert to pseudo_coords to dataframe
  bias_points <- coords2df(bias_points)
  
  # add weights column to pseudo_coords
  bias_points$weight <- rep(weight, nrow(bias_points))
  
  return(bias_points)
  
}


# remove NA values from a ppm object
ppmClean <- function (ppm) {

  # remove any points that are NA and issue a warning
  if (any(is.na(ppm))) {

    # copy the original data
    ppm_old <- ppm
    # remove NAs
    ppm <- na.omit(ppm)

    # report which ones were removed
    rm <- as.vector(attributes(ppm)$na.action)
    which_rm <- ppm_old$points[rm]

    if (any(which_rm == 0)) {
      warning(sprintf('Removed %i integration points for which covariate values could not be assigned.
                      This may affect the results of the model, please check alignment between covariates and area.',
                      sum(which_rm == 0)))
    }

    if (any(which_rm == 1)) {
      warning(sprintf('Removed %i observed points for which covariate values could not be assigned.
                      This may affect the results of the model, please check alignment between covariates and coords.',
                      sum(which_rm == 1)))
    }

  }

  return (ppm)

}
