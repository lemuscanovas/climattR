library(sf)
library(terra)
library(tidyverse)

as_ts <- function(x, detrend = F,hour2day_fun = "mean", domain = NULL, sf_obj = NULL) {
  
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
    # linear detrend
    message("Starting to detrend the full time series")
    yr <- as.factor(year(time_dat))
    dat_yr_mean <- tapp(dat,yr,"mean")
    dat_yr_mean[is.na(dat_yr_mean)] <- -9999 #avoiding problems with NAs out of region
    
    lyrs <- 1:nlyr(dat_yr_mean)
    
    .remove_slope <- function(x) { 
      m <- lm(x ~ lyrs) 
      rm_slope <- residuals(m) + predict(m)[1] # elimina slope
      return(rm_slope)
    }
    
    # remove slope from yearly time series
    dat_yr_rm_slope <- app(dat_yr_mean, .remove_slope) 
    dat_yr_mean[dat_yr_mean <  -9998] <- NA
    dat_yr_rm_slope[dat_yr_rm_slope < -9998] <- NA
    # getting slope for each year
    slope_yr <- dat_yr_mean - dat_yr_rm_slope
    slope_yr[[1]][slope_yr[[1]]< Inf] <- 0 # first no changes
    # remove slope from daily ts
    yrs <- unique(year(time_dat))
    ind <- seq_along(yrs)
    
    .detrend_fun <-  function(ii){
      ind_yr <- which(year(time_dat)==yrs[ii])
      selyr_dat <- dat[[ind_yr]]
      selyr_dat_slope <- slope_yr[[ii]]
      selyr_dat_detrended <- selyr_dat - selyr_dat_slope
      
      return(selyr_dat_detrended)
    }
    
    # detrended
    dat <- pblapply(ind, .detrend_fun) %>% rast()
  }
  
  if(!is.null(domain)){
    domain <- ext(domain) %>% as.polygons()
  o <- terra::extract(dat,domain, fun = mean, na.rm = T ) %>%
    t() %>% as.vector()
  }
  if(!is.null(sf_obj)){
    o <- terra::extract(dat,vect(sf_obj), fun = mean, na.rm = T ) %>%
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
