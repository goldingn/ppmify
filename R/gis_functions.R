# convert an area object into a raster with each cell giving the area if area is
# either rasterize or resample to make sure it has the same layout as the
# template
areaRaster <- function (area, template) {

  # convert area into a raster if needed
  if (inherits(area, c('SpatialPolygons', 'SpatialPolygonsDataFrame'))) {
    area <- rasterize(area, template)
  } else {
    # otherwise resample
    area <- resample(area * 0, template, method = 'ngb')
  }

  # get cell areas
  if (isLonLat(area)) {
    # cell-wise areas in km ^ 2
    cell_area <- area(area)
  } else {
    # otherwise, assume cells are equal area
    # get overall surface area in square kilometres
    surface_area <- prod(dim(area)[1:2] * res(area)) * 0.1 ^ 6
    # set area in each cell
    cell_area <- setValues(area, surface_area / ncell(area))

  }

  # mask out areas outside the area of interest
  cell_area <- mask(cell_area, area)

  return (cell_area)

}

# given a grid raster and an area object, calculate the surface area of area
# falling in each cell of grid. res sets the resolution of the area estimate,
# relative to the grid
gridArea <- function (grid, area, res = 10) {

  # disaggregate grid 10x
  grid_disagg <- disaggregate(grid, res)
  projection(grid_disagg) <- projection(grid)

  # convert area into a raster of the same size, giving cell areas
  area_ras <- areaRaster(area, grid_disagg)

  # sum by cell number
  sums <- zonal(area_ras, grid_disagg, fun = 'sum')

  # put back into grid
  grid_area <- reclassify(grid, sums)

  return (grid_area)

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

# calculate the dimensions of an extent in km (either an extent object or
# four-element vector in the right order), either in projected or spherical
# space
extentDim <- function (extent, lonlat =  TRUE) {
  # coerce to vector if necessary
  if (inherits(extent, 'extent')) extent <- as.vector(extent)
  if (lonlat) {
    # dimensions in degrees
    height <- abs(diff(extent[1:2]))
    width <-  abs(diff(cos(extent[3:4])))
    # Scaling to get spherical surface area in km2
    scaling <- (6371 ^ 2 * pi) / 180
    surface_area <- width * height * scaling
    # ratio between GCD height and width
    ratio <- lonLatRatio(extent)
    # calculate equivalent dimensions in km
    w <- sqrt(surface_area / ratio)
    dim <- c(w, w * ratio)
  } else {
    # else assume a rectangle in m and convert to km
    dim <- abs(diff(extent)[c(1, 3)]) * 0.1 ^ 3
  }
  return (dim)
}

# get ratio between height and width in great circle distance, given an extent
# vector in lat/long
lonLatRatio <- function (extent) {
  # lower left point
  p1 <- matrix(extent[c(1, 3)], nrow = 1)
  # upper left and lower right points
  p2 <- rbind(extent[c(1, 4)], extent[c(2, 3)])
  # get ratio between distances
  dists <- pointDistance(p1, p2, lonlat = TRUE)
  ratio <- dists[1] / dists[2]
  return (ratio)
}
