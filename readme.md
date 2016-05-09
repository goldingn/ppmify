# ppmify

### convert your Poisson regression into a Poisson point process model

[![Build Status](https://travis-ci.org/goldingn/ppmify.svg)](https://travis-ci.org/goldingn/ppmify)
[![codecov.io](https://codecov.io/github/goldingn/ppmify/coverage.svg?branch=master)](https://codecov.io/github/goldingn/ppmify)
[![cran version](http://www.r-pkg.org/badges/version/ppmify)](http://cran.rstudio.com/web/packages/ppmify)

**package under construction - not ready to use just yet!**

![construction](https://camo.githubusercontent.com/4a7cf94aedbd23c13cc2d75fdc3b2af5c816c208/687474703a2f2f7374617469632e646967672e636f6d2f7374617469632f696d616765732f6469676765722e676966)

**ppmify** (pronounced *p-p-m-if-eye*) is a micropackage to help set up Poisson point process models (PPMs) for point data.
PPMs can be fitted using standard software for Poisson regression, after a little bit of fiddling with the data (adding integration points as pseudo-observations and calculating integration weights).
This process isn't too difficult, but can be a bit time consuming.
The function `ppmify()` takes the data needed for fitting a PPM, does this fiddling for you, and then returns a dataframe-like object that you can use in your favourite Poisson modelling software.
That could be GLM, a generalised boosted model (AKA boosted regression tree), elasticnet regression, Gaussian process regression or whatever.

ppmify aims to make the data preparation step a bit easier. 
At some point in the future it should be able to set up slightly more complex integration methods, thinned PPMs, PPMs with attraction/repulsion between points, and possibly even multi-species PPMS - watch this space!

However, this isn't an end-to end solution.
Once `ppmify()` has handed over the processed data, it's up to the user to fit the Poisson model in the right way.
We'll provide some guidance here, but would advise users to familiarise themselves with PPMs and this modelling approach first. 

