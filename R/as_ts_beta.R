library(sf)
library(terra)
library(tidyverse)

as_ts <- function(x, detrend = F, k = 2,hour2day_fun = "mean", domain = NULL, sf_obj = NULL) {
  
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
  if(length((time(dat))) != length(unique(time(dat)))){
    if(hour2day_fun != "mean"){
      hour2day_FUN <- match.fun(FUN = hour2day_fun)
      dat <- dat %>%
        tapp(dat, as.factor(time_dat),hour2day_FUN, na.rm = T)
    }
    dat <- tapp(dat, as.factor(time_dat),"mean", na.rm = T)
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
    domain <- ext(domain) %>% as.polygons()
  o <- terra::extract(dat,domain, fun = mean, na.rm = T ) %>%
    t() %>% as.vector()
  }
  if(!is.null(sf_obj)){
    o <- terra::extract(dat,vect(sf_obj), fun = NULL, na.rm = T,method = "simple") %>%
    t() %>% as.vector()
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
