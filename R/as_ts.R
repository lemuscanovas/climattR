#' Convert Spatial Raster Data to Time Series
#'
#' This function processes a spatial raster dataset, optionally detrending and aggregating by time,
#' to convert it into a time series format. It supports conversion from hourly to daily data and can
#' apply a domain or spatial filter using an `sf` object.
#'
#' @param x A file path to a NetCDF file or a `SpatRaster` object.
#' @param detrend Logical; should the data be detrended? Defaults to FALSE.
#' @param k Degree of the polynomial for detrending if `detrend` is TRUE.
#' @param hour2day_fun A character string specifying the function to aggregate hourly data to daily. Default is "mean".
#' @param domain Optional; a spatial object or extent to restrict the analysis within a specific domain.
#' @param sf_obj Optional; an `sf` object used to extract data points based on spatial features.
#'
#' @return A `tibble` containing time series data with two columns: time and the variable of interest.
#'
#' @importFrom terra rast time app extract
#' @importFrom lubridate as_date
#' @importFrom tibble tibble 
#' @importFrom pracma polyfit polyval
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#' TX <- terra::rast(system.file("extdata", "tx_0608_1950_2023.nc", package = "climattR"))
#' ts <- as_ts(TX, detrend = FALSE, hour2day_fun = "mean")}
#' @export

as_ts <- function(x, 
                  detrend = FALSE, 
                  k = 2, 
                  hour2day_fun = "mean", 
                  domain = NULL, 
                  sf_obj = NULL) {
  
  if(class(x)[1] == "SpatRaster"){
    dat <- x
  }else{
    dat <- rast(x)  
  }
  
  time_dat <- as_date(time(dat))
  
  if(is.Date(time_dat) == F){
    stop("time is not a date object. Please, check time values in your data")
  }
  
  # If hourly, converts to daily
  if(length(time(dat)) != length(unique(time(dat)))){
    dat <- tapp(dat, as.factor(time_dat), hour2day_fun, na.rm = TRUE)
  }    
  
  
  if(isTRUE(detrend)){
    detrend_pracma <- function(y, k) {
      fit <- polyfit(seq_along(y), y, k)
      y_pred <- polyval(fit, seq_along(y))
      y - y_pred
    }
    
    dat <- app(dat, detrend_pracma, k = k)
  }
  
  
  if(!is.null(domain)){
    domain_poly <- ext(domain) %>% as.polygons()
    o <- terra::extract(dat, domain_poly, fun = mean, na.rm = TRUE) %>%
      t() %>% as.vector()
  } else if(!is.null(sf_obj)){
    o <- terra::extract(dat, vect(sf_obj), fun = NULL, na.rm = TRUE, method = "simple") %>%
      t() %>% as.vector()
  } else {
    # No spatial filter: spatial mean over all cells
    o <- c(NA_real_, terra::global(dat, "mean", na.rm = TRUE)[[1]])
  }
  
  dates <- as_date(time(dat))
  
  if(is.Date(dates) == F){
    stop("time is not a date object. Please, check time values in your data")
  }
  
  ts <- tibble(dates, o[-1]) %>% 
    setNames(c("time", "var")) %>%
    group_by(time) %>%
    summarise(var = mean(var, na.rm = T),.groups = "drop") %>%
    setNames(c("time", varnames(dat)))
  
  
  return(ts) 
}
