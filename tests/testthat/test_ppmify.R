context('ppmify')


test_that('coords handles expected inputs', {

  # input datatypes
  mx <- matrix(cbind(runif(10),
                     runif(10)),
               ncol = 2)

  df <- data.frame(mx)
  sp <- SpatialPoints(mx)
  spdf <- SpatialPointsDataFrame(mx, df)

  # run models
  ppm_matrix <- ppmify(coords = mx)
  ppm_df <- ppmify(coords = df)
  ppm_sp <- ppmify(coords = sp)
  ppm_spdf <- ppmify(coords = spdf)

  # check outputs
  outclass <- c('data.frame', 'ppm')
  expect_is(ppm_matrix, outclass)
  expect_is(ppm_df, outclass)
  expect_is(ppm_sp, outclass)
  expect_is(ppm_spdf, outclass)

})

test_that('coords handles unexpected inputs', {

  # some input datatypes
  ch <- 'foo'
  int <- 1L
  num <- 3.1
  f <- a ~ b

  # force errors
  expect_error(ppmify(coords = ch))
  expect_error(ppmify(coords = int))
  expect_error(ppmify(coords = num))
  expect_error(ppmify(coords = f))

})
