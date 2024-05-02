climattR <img src="img/logo.png" align="right" alt="" width="140" />
=========================================================
# `climattR`: Rapid climate extreme event attribution for regional and local areas


[![CRAN status](https://www.r-pkg.org/badges/version/climattR)](https://cran.r-project.org/package=climattR)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

## Overview

**climattR** is a toolbox for performing rapid attribution of extreme weather events using the flow analogs approach. *This package is still under developement*

----

Contents:

* [Why this package](#why-this-package)
* [How it works](#how-it-works)
  * [Installation](#installation)
* [Package citation](#package-citation)
* [Contact](#contact)

----

## Why this package
To date there is no R package available which have implemented the simple but effective flow analogues methodology to perform extreme weather event attribution. `climattR` contains a set of functions to 1) prepare data to employ flow analogues; 2) Find flow analogue days to a target date; 3) Reconstruct surface fields for extreme weather events in past and present (also future) climate conditions;

## How it works

### Installation

``` r
# To install the latest version from Github:
# install.packages("remotes")
# remotes::install_github("lemuscanovas/climattR")

library(climattR)
```

### First steps. Rapid attribution of 

``` r
z <- list.files(path = "inst/testdata/",pattern = "geopot",full.names = T) %>% rast()
t2m <- list.files(path = "inst/testdata/",pattern = "2m_t",full.names = T) %>% rast()
```

``` r
## Preparing data for analogs
dates <- seq(as_date("2022-06-12"), as_date("2022-06-21"),"day")

dat4an <- prepare_data(x = z,level = NULL,
                      event_dates =  dates,time_window = 31)

```


## Package citation

Using synoptReg for research publication?  Please **cite it**! I'm an early career scientist and every citation matters.

***Lemus-Canovas, M., Lopez-Bustins, J.A., Martin-Vide, J., Royé, D.***, 2019. *synoptReg: An R package for computing a synoptic climate classification and a spatial regionalization of environmental data*. Environmental Modelling & Software, Vol. 118,114-119pp, ISSN 1364-8152, https://doi.org/10.1016/j.envsoft.2019.04.006

## Contact

Feel free to contact me: mlemus@ub.edu
