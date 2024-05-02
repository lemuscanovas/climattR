#' Prepare data from NetCDF with spatial and temporal processing
#'
#' This function processes NetCDF data with options for selecting specific pressure levels,
#' detrending, scaling, and extracting specific time windows around provided event dates.
#' It supports handling both daily and sub-daily resolutions.
#'
#' @param x A file path to a NetCDF file or a `SpatRaster` object.
#' @param level Optional; specify a particular pressure level to extract.
#' @param event_dates Dates of specific events to focus on within the data.
#' @param time_window Number of days before and after event_dates to include in the output.
#' @param analog_months Optional; specify months to be used in analysis if no specific time window provided.
#' @param detrend Logical; should the data be detrended? Defaults to FALSE.
#' @param k Degree of the polynomial for detrending if `detrend` is TRUE.
#' @param scale Logical; should the data be scaled? Defaults to FALSE.
#'
#' @return A list with two elements: `ts_wo_event`, the timeseries data excluding the event dates and window,
#'         and `event`, the timeseries data for the event dates.
#'
#' @importFrom terra rast time app
#' @importFrom lubridate as_date year month first last
#' @importFrom tidyverse %>% tibble filter
#' @importFrom pbapply pbsapply
#' @importFrom stringr str_detect str_sub str_c
#' @importFrom pracma polyfit polyval
#' @importFrom stats scale
#' @examples
#' nc_file <- system.file("extdata", "example.nc", package = "yourPackageName")
#' event_dates <- as.Date(c("2022-06-01"))
#' result <- prepare_data(nc_file, event_dates = event_dates, time_window = 30)
#' @export

prepare_data <- function(x, level = NULL, event_dates,
                         time_window = 31, analog_months = NULL,
                         detrend = F, k = 2, scale = F) {
  # function implementation remains the same

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
  
  # seq_event <- seq(min(as_date(event_dates))-time_window, 
  #          max(as_date(event_dates))+time_window, 
  #          by = "day")
  
  if(!is.null(time_window)){
  .time_windows = function(x){
    seq(min(as_date(event_dates))-time_window, 
        max(as_date(event_dates))+time_window, 
        by = "day") %>%
      str_sub(start = -5,end = -1) %>%
      str_c(x,sep = "-") %>% 
      as.vector()
    }
    time_window_an <-  sapply(yr_seq, FUN = .time_windows) %>% mdy()

  }
    
    if(is.null(time_window) & !is.null(analog_months )){
      
      dates <- seq(as_date(paste0(year_range[1],"-01-01")), 
          as_date(paste0(year_range[2],"-12-31")), 
          by = "day")
      mo = which(month(dates) %in% analog_months)
      
      time_window_an <- dates[mo]
    }
  
  
  message("Additionally, if hourly data  provided then it is converted to daily mean.")
  yr_seq <-  yr_seq[!yr_seq %in% event_yr]
  
  time_all <- filter(ts_nc, time %in% time_window_an)
  nc_timeseries_dm_an <- nc_dayhour[[time_all$id]] 
  
  nc_timeseries_dm_an <- nc_timeseries_dm_an %>%
    tapp(as.factor(time_all$time),"mean")
  
  time_dy <- as_date(time_all$time) %>% unique()
  time(nc_timeseries_dm_an) <- time_dy
 
  if(isTRUE(detrend)){
    detrend_pracma <- function(y, k) {
      fit <- polyfit(seq_along(y), y, k)
      y_pred <- polyval(fit, seq_along(y))
      y - y_pred
    }
    
    nc_timeseries_dm_an <- app(nc_timeseries_dm_an, detrend_pracma,k=k)
    # dat <- c(dat[[1]],dat)
  }
  
  if(isTRUE(scale)){
    nc_timeseries_dm_an <- app(nc_timeseries_dm_an,scale.default)
  }
  
  time4analogs <- unique(time_all$time)
  nc_event_dm_an <- nc_timeseries_dm_an[[which(time4analogs %in% as_date(event_dates))]]
  
  # Excluding event dates and window dates of the event
  
  if(!is.null(time_window)){
  windows_event <- seq(min(as_date(event_dates))-time_window, 
                       max(as_date(event_dates))+time_window, 
                       by = "day")
  }
  
  if(is.null(time_window) & !is.null(analog_months )){
    
    dates <- seq(as_date(paste0(year_range[1],"-01-01")), 
                 as_date(paste0(year_range[2],"-12-31")), 
                 by = "day")
    mo = which(month(dates) %in% analog_months)
    monthly_dates <- dates[mo]
    
    # tenim un episodi que abraça dos anys?
    dos_anys <- ifelse(analog_months %in% c(12,1), T,F)
    if(sum(dos_anys)> 1){
      event_yr <-  c(event_yr-1,event_yr)
      windows_event <- seq.Date(as_date(paste0(event_yr[1],"-",analog_months[1],"-01")), 
                                as_date(paste0(event_yr[2],"-",analog_months[length(analog_months)],"-31")), "day")
    } else{
      windows_event <- seq.Date(as_date(paste0(event_yr,"-",analog_months[1],"-01")),
                                as_date(paste0(event_yr,"-",analog_months[length(analog_months)],"-31")), "day")
    }
    
  }
  # windows_event <- seq(min(as_date("2022-06-01")), 
  #                      max(as_date("2022-08-31")), 
  #                      by = "day")
  
  if(detrend == T){
    nc_timeseries_dm_an <- nc_timeseries_dm_an[[which(!time4analogs %in% windows_event)]]
    time_all <- time4analogs[which(!time4analogs %in% windows_event)]
    
  } else{
  
  nc_timeseries_dm_an <- nc_timeseries_dm_an[[which(!time4analogs %in% windows_event)]]
  time_all <- time4analogs[which(!time4analogs %in% windows_event)]
  }
  
  terra::time(nc_timeseries_dm_an) <- as.POSIXlt(as_date(time_all))
  names(nc_timeseries_dm_an) <- as.POSIXlt(as_date(time_all))
  
  terra::time(nc_event_dm_an) <- as.POSIXlt(as_date(event_dates))
  names(nc_event_dm_an) <- as.POSIXlt(as_date(event_dates))
  
  
return(list(ts_wo_event = nc_timeseries_dm_an, 
            event = nc_event_dm_an))
}
