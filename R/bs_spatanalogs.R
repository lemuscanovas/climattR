#' Bootstrapping Spatial Analogs
#'
#' This function performs bootstrapping on spatial analogs for a given set of event dates.
#' It allows for detrending, anomaly calculation, and custom event functions on a `SpatRaster`.
#' The function returns a set of spatial datasets simulating possible scenarios.
#'
#' @param x A file path to a NetCDF file or a `SpatRaster` object.
#' @param analogs A data frame containing analog dates and corresponding period labels.
#' @param n The number of bootstrap samples to generate.
#' @param event_fun A character string specifying the aggregation function to apply, default is "mean".
#' @param anom Logical; if TRUE, calculates anomalies based on a reference period.
#' @param ref_period A numeric vector of length two specifying the start and end years for the reference period, required if `anom` is TRUE.
#' @param replace Logical; should sampling be with replacement? Defaults to TRUE.
#' @param detrend Logical; should the data be detrended? Defaults to FALSE.
#' @param k Degree of the polynomial for detrending if `detrend` is TRUE.
#'
#' @return A list of `SpatRaster` objects for each unique period in the analogs, representing simulated bootstrap samples.
#'
#' @importFrom terra rast app time
#' @importFrom lubridate as_date year
#' @importFrom pracma polyfit polyval
#' @importFrom magrittr %>% 
#' @importFrom dplyr inner_join filter slice_sample group_by mutate select
#' @importFrom pbapply pblapply
#' @examples
#' dat_path <- system.file("extdata", "example.nc", package = "yourPackageName")
#' analogs_data <- data.frame(time = as.Date('2021-01-01') + 0:10, period = rep('2021', 11))
#' result <- bs_spatanalogs(dat_path, analogs_data, n = 100, event_fun = "mean", anom = TRUE, ref_period = c(1980, 2000))
#' @export

bs_spatanalogs <- function(x, analogs, n = 1000, 
                           event_fun = "mean",
                           anom = NULL, 
                           ref_period = NULL,
                           replace = TRUE,
                           detrend = FALSE, k = 2) {
  
  event_FUN <- match.fun(FUN = event_fun)
  
  # Initialize raster
  nc_dayhour <- if (inherits(x, "SpatRaster")) {
    dat <- x
  } else {
    dat <- rast(x)
  }
  
  dat <- dat * 1 # stored in disk (?)
  
  time_dat <- as_date(time(dat)) # if hourly it will be converted to daily
  
  if(is.Date(time_dat) == F){
    stop("time is not a date object. Please, check time values in your data")
  }
  
  # If hourly, converts to daily
  message("Caution! if hourly data is provided then it is converted into daily mean.")
  
  if(length((time(dat))) != length(unique(time(dat)))){
      dat <- dat %>% tapp(dat, as.factor(time_dat),"mean", na.rm = T)
    }
  
  # Apply detrending if specified
  if(isTRUE(detrend)){
    detrend_pracma <- function(y, k) {
      fit <- polyfit(seq_along(y), y, k)
      y_pred <- polyval(fit, seq_along(y))
      y - y_pred
    }
    
    dat <- app(dat, detrend_pracma, k = k)
  }
  
  if(isTRUE(anom) & length(ref_period == 2)){
    
    ref_mean <- dat %>% 
      subset(which(year(time_dat) >= ref_period[1] & year(time_dat) <= ref_period[2])) %>%
      app("mean", na.rm =T)
    
    
  }else if(isFALSE(anom)){
    
  } else if(isTRUE(anom) & is.null(ref_period)){
    stop("Please, a reference period is required. eg: ref_period = c(1981,2010)")
    
  }else{
    stop("Please, the reference period is as follows: ref_period = c(1981,2010)")
  }

  ts_nc <-tibble(id = seq_along(time(dat)), # id to subset
                 time = time_dat) 
  
  # Calculation of bootstrapped sd and mean for the different periods set ------
  yr_split <- analogs$period %>% 
    unique() 
  
  sim_bs_l <- list()
  for(ii in seq_along(yr_split)){
    
    # selecting  analogs
    dat_subset <- inner_join(analogs, ts_nc, by = "time") %>%
    filter(str_detect(period, yr_split[ii]))
    
      bootstrap <- function(x){
      cases <- dat_subset %>% group_by(time_obj) %>%
        slice_sample(n = 1, replace = replace)
    
      # Computing mean
      sim_event <- dat[[cases$id]] %>% 
        app(event_fun) # summarising the event
  
    }

message(paste0("bootstrapping the event ", n," times in the ", yr_split[ii], ":"))
      
  system.time(sim_bs<- pbapply::pblapply(1:n, bootstrap, cl = 1) %>% 
                rast())
  names(sim_bs) <- paste0("sim",1:n)
  
  if(isTRUE(anom)){
    sim_bs <- sim_bs - project(ref_mean, sim_bs)
  }
  
  sim_bs_l[[ii]] <- sim_bs
  }

  sim_bs_sds <- sds(sim_bs_l)
  names(sim_bs_sds) <- yr_split
  
  return(sim_bs_sds)

  }
