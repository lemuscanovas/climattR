library(tidyverse)
library(lubridate)
library(terra)

spat_analogs <- function(x, analogs_subperiods,detrend = F, 
                         conversion_fun = NULL, hour2day_fun = "mean" ,event_fun, save = F){
  
  event_FUN <- match.fun(FUN = event_fun)

  
  # Reading analogs dates and VOI nc ----------------------------------------
  if(class(x)[1] == "SpatRaster"){
    dat <- x
  }else{
    dat <- rast(x)  
  }
  
  varname <- varnames(dat)
  
  time_dat <- as_date(time(dat)) # if hourly it will be converted to daily
  
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
    dat_yr_mean <- tapp(dat,yr,hour2day_fun)
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
  
  if(!is.null(conversion_fun)){
    conversion_FUN <- match.fun(FUN = conversion_fun)
    dat <- dat %>%
      app(conversion_FUN)
  }
  

  ts_nc <-tibble(id = seq_along(time(dat)), # id to subset
                 time = time_dat) 
  
  # Computing sd and mean for a counterfactual world -------------------------------
  
  # to select first period (counterfactual)
  yr_split <- analogs_subperiods$period %>% 
    str_sub(start = -4,end = -1) %>% 
    as.numeric() %>% 
    min %>% 
    as.character()
  
  # selecting factual analogs
  counter <- inner_join(analogs_subperiods, ts_nc, by = "time") %>%
    filter(str_detect(period, yr_split))
  
  # Computing mean
  counter_daybyday_mean <- dat[[counter$id]] %>% # selectinc VOI analogs
    tapp(as.factor(counter$time_obj),"mean")#mean

  # if the event is greater than 1 day, then compute mean or sum for the whole event
    if(length(unique(analogs_subperiods$time_obj))>1){
      counter_event_mean <- app(counter_daybyday_mean,event_FUN) # sum days pcp (other fun for other vars)
    }
  
  # computing sd
  counter_daybyday_sd <- dat[[counter$id]] %>% 
    tapp(as.factor(counter$time_obj),"sd")

  if(length(unique(analogs_subperiods$time_obj))>1){
    # mean days sd pcp (other fun for other vars)
    counter_event_sd <- app(counter_daybyday_sd,"mean") 
    }
  
  # Computing sd and mean for a factual world -------------------------------
  
  # to select second period (factual)
  yr_split <- analogs_subperiods$period %>% 
    str_sub(start = -4,end = -1) %>% 
    as.numeric() %>% 
    min %>% sum(1) %>%
    as.character()
  
  # selecting counterfactual analogs
  factual <- inner_join(analogs_subperiods, ts_nc, by = "time") %>%
    filter(str_detect(period, yr_split))
  
  # Computing mean
  factual_daybyday_mean <- dat[[factual$id]] %>% 
    tapp(as.factor(factual$time_obj),"mean")

    # if the event is greater than 1 day, then compute mean or sum for the whole event
    if(length(unique(analogs_subperiods$time_obj))>1){
    
      # sum days pcp (other fun for other vars)
      factual_event_mean <- app(factual_daybyday_mean,event_FUN) 
      
    }
  
  factual_daybyday_sd<- dat[[factual$id]] %>% 
    tapp(as.factor(factual$time_obj),"sd")

    if(length(unique(analogs_subperiods$time_obj))>1){
    # sum days pcp (other fun for other vars)
      factual_event_sd <- app(factual_daybyday_sd,"mean") 
    }
  
  
  # Saving objects ----------------------------------------------------------
    if(isTRUE(save)){
      folder_name <- "01_counter_factual_spatanalogs"
      x <- dir.create(file.path(folder_name), showWarnings = FALSE)      
      terra::time(factual_daybyday_mean) <- as.POSIXlt(unique(analogs_subperiods$time_obj))
      writeCDF(factual_daybyday_mean,
               filename = paste0(folder_name,"/", varname,"_factual_daily_mean.nc"),
               varname = varnames(dat), zname = "time", 
               prec = "float", 
               overwrite = T,compression = 3)
      
      terra::time(factual_daybyday_sd) <- as.POSIXlt(unique(analogs_subperiods$time_obj))
      writeCDF(factual_daybyday_sd,
               filename = paste0(folder_name,"/",varname,"_factual_daily_sd.nc"),
               varname = varnames(dat), zname = "time", 
               prec = "float", 
               overwrite = T,compression = 3)
      
      terra::time(counter_daybyday_mean) <- as.POSIXlt(unique(analogs_subperiods$time_obj))
      writeCDF(counter_daybyday_mean,
               filename = paste0(folder_name,"/",varname,"_counterfactual_daily_mean.nc"),
               varname = varnames(dat), zname = "time", 
               prec = "float", 
               overwrite = T,compression = 3)
      
      terra::time(counter_daybyday_sd) <- as.POSIXlt(unique(analogs_subperiods$time_obj))
      writeCDF(counter_daybyday_sd,
               filename = paste0(folder_name,"/",varname,"_counterfactual_daily_sd.nc"),
               varname = varnames(dat), zname = "time", 
               prec = "float", 
               overwrite = T,compression = 3)
      if(length(unique(analogs_subperiods$time_obj))>1){
        
        writeCDF(factual_event_mean,
                 filename = paste0(folder_name,"/",varname,"_factual_event_mean.nc"),
                 varname = varnames(dat), zname = "time", 
                 prec = "float", 
                 overwrite = T,compression = 3)
        
        writeCDF(counter_event_mean,
                 filename = paste0(folder_name,"/",varname,"_counterfactual_event_mean.nc"),
                 varname = varnames(dat), zname = "time", 
                 prec = "float", 
                 overwrite = T,compression = 3)
        
      }
      
    }
  
  if(length(unique(analogs_subperiods$time_obj))>1){
    return(list(factual_all_mean = factual_daybyday_mean, 
                factual_all_sd = factual_daybyday_sd,
                counter_all_mean = counter_daybyday_mean,
                counter_all_sd = counter_daybyday_sd,
                factual_event_mean = factual_event_mean,
                counter_event_mean = counter_event_mean))
    }else{
      return(list(factual_event_mean = factual_daybyday_mean,
                  factual_event_sd = factual_daybyday_sd,
                  counter_event_mean = counter_daybyday_mean,
                  counter_event_sd = counter_daybyday_sd))
    }
}

