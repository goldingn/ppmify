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

  # convert area into a raster of the same size
  cell_area <- areaRaster(area, grid_disagg)

  # sum by cell number
  sums <- zonal(cell_area, grid_disagg, fun = 'sum')

  # put back into grid
  grid_area <- reclassify(grid, sums)

  return (grid_area)

}
