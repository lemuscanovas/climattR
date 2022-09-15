library(tidyverse)
library(lubridate)
library(terra)
library(pbapply)


tidync_attr <- function(x, level = NULL, detrend = F, scale =  F, extent = NULL, aggregate = NULL, 
                        rotate = F,event_dates, time_window, save){

  # Reading analogs dates and VOI nc ----------------------------------------
  if(class(x)[1] == "SpatRaster"){
    nc_dayhour <- x
  }else{
    nc_dayhour <- rast(x)  
  }
  
  if(is.numeric(level) == T){ # selecting pressure level
    
  levs <- names(nc_dayhour)
  l <- which(str_detect(levs, paste0(level,"_")))
  
  nc_dayhour <- nc_dayhour[[l]] 
  }
  
  ts_daily  = as_date(terra::time(nc_dayhour)) %>% sort()
  ts_yearly = year(ts_daily) %>% unique()
  year_range <- c(first(ts_yearly),last(ts_yearly))
  
  if(!is.Date(ts_daily)){
    stop("NetCDF provided has uncodified dates. Check time variable")
  } 
  
  ts_nc <- tibble(id = seq_along(ts_daily),
                      time = ts_daily)
  
  event_yr = year(event_dates[1])
  yr_seq <- seq(year_range[1],year_range[2]) 
  
  seq_event <- seq(min(as_date(event_dates))-time_window, 
           max(as_date(event_dates))+time_window, 
           by = "day")
  
  .time_windows = function(x){
    seq(min(as_date(event_dates))-time_window, 
        max(as_date(event_dates))+time_window, 
        by = "day") %>%
      str_sub(start = -5,end = -1) %>%
      str_c(x,sep = "-") %>% 
      as.vector()
  }
  
  time_window_an <-  sapply(yr_seq, FUN = .time_windows) %>% mdy()
  
  message("Resampling and cropping full time series. Hourly data converted to daily.")
  message("It can take a while...")
  yr_seq <-  yr_seq[!yr_seq %in% event_yr]
  
  time_all <- filter(ts_nc, time %in% time_window_an)
  nc_timeseries_dm_an <- nc_dayhour[[time_all$id]] 
  
  if(isTRUE(rotate)){
    nc_timeseries_dm_an <- terra::rotate(nc_timeseries_dm_an)
  } 
  
  if(is.null(extent)){
    extent <- ext(nc_timeseries_dm_an)
  }
  
  
  if(is.null(aggregate)){
    nc_timeseries_dm_an <- nc_timeseries_dm_an %>%
      crop(extent) %>% 
    tapp(as.factor(time_all$time),"mean")
  } else{
    nc_timeseries_dm_an <- nc_timeseries_dm_an %>%
      aggregate(aggregate,"mean") %>%
      crop(extent) %>% 
      tapp(as.factor(time_all$time),"mean")
  }
  
  time_dy <- as_date(time_all$time) %>% unique()
  if(isTRUE(detrend)){
    # linear detrend
    message("Starting to detrend the full time series")
    yr <- as.factor(year(time_dy))
    nc_timeseries_dm_an_yr <- tapp(nc_timeseries_dm_an,yr,"mean")
     
    lyrs <- 1:nlyr(nc_timeseries_dm_an_yr)
    
    .remove_slope <- function(x) { 
      m <- lm(x ~ lyrs) 
      rm_slope <- residuals(m) + predict(m)[1] # elimina slope
      return(rm_slope)
    }
    
    # remove slope from yearly time series
    nc_timeseries_dm_an_yr_rm_slope <- app(nc_timeseries_dm_an_yr, .remove_slope) 
    
    # getting slope for each year
    slope_yr <- nc_timeseries_dm_an_yr - nc_timeseries_dm_an_yr_rm_slope
    
    # remove slope from daily ts
    yrs <- unique(year(time_dy))
    ind <- seq_along(yrs)
    
    .detrend_fun <-  function(ii){
      ind_yr <- which(year(time_dy)==yrs[ii])
      selyr_ncfull_ts <- nc_timeseries_dm_an[[ind_yr]]
      selyr_ncfull_slope <- slope_yr[[ii]]
      selyr_ncfull_detrended <- selyr_ncfull_ts - selyr_ncfull_slope
      
      return(selyr_ncfull_detrended)
    }
    
    # detrended
    nc_timeseries_dm_an <- pblapply(ind, .detrend_fun) %>% rast()
  }
  
  if(isTRUE(scale)){
    nc_timeseries_dm_an <- scale(nc_timeseries_dm_an)
  }
  
  time4analogs <- unique(time_all$time)
  nc_event_dm_an <- nc_timeseries_dm_an[[which(time4analogs %in% as_date(event_dates))]]
  
  # Excluding event dates and window dates of the event
  windows_event <- seq(min(as_date(event_dates))-time_window, 
                       max(as_date(event_dates))+time_window, 
                       by = "day")
  
  nc_timeseries_dm_an <- nc_timeseries_dm_an[[which(!time4analogs %in% windows_event)]]
  
  time_all <- time4analogs[which(!time4analogs %in% windows_event)]
  terra::time(nc_timeseries_dm_an) <- as.POSIXlt(unique(time_all))
  terra::time(nc_event_dm_an) <- as.POSIXlt(as_date(event_dates))
  
  if(isTRUE(save)){
    folder_name <- "00_data4analogs"
    x <- dir.create(file.path(folder_name), showWarnings = FALSE)
    writeCDF(nc_timeseries_dm_an,
             filename = paste0(folder_name,"/",varnames(nc_event_dm_an),level,"_",
                               paste0(year_range,collapse = ""),
                               "_woevent_day.nc"),
             varname = varnames(nc_event_dm_an), zname = "time", 
             prec = "float", 
             overwrite = T,compression = 3)
    
    writeCDF(nc_event_dm_an,
             filename = paste0(folder_name,"/",varnames(nc_event_dm_an),level,"_",
                               paste0(year_range,collapse = ""),
                               "_event_day.nc"),
             varname = varnames(nc_event_dm_an), zname = "time", 
             prec = "float", 
             overwrite = T, compression = 3)
  }
  
return(list(ts_wo_event = nc_timeseries_dm_an, 
            event = nc_event_dm_an))
}
