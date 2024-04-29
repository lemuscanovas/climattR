library(tidyverse)
library(lubridate)
library(pracma)

bs_spatanalogs<- function(x,analogs,n = 1000, 
                          event_fun = "mean",
                          anom = NULL, 
                          ref_period = NULL,
                          replace = T,
                          detrend = F,k = 2){
  
  event_FUN <- match.fun(FUN = event_fun)
  
  # Reading analogs dates and VOI nc ----------------------------------------
  if(class(x)[1] == "SpatRaster"){
    dat <- x
  }else{
    dat <- rast(x)  
  }
  
  time_dat <- as_date(time(dat)) # if hourly it will be converted to daily
  
  if(is.Date(time_dat) == F){
    stop("time is not a date object. Please, check time values in your data")
  }
  
  # If hourly, converts to daily
  message("Caution! if hourly data is provided then it is converted into daily mean.")
  
  if(length((time(dat))) != length(unique(time(dat)))){
      dat <- dat %>%
        tapp(dat, as.factor(time_dat),"mean", na.rm = T)
    }
  
  if(isTRUE(detrend)){
    detrend_pracma <- function(y, k) {
      fit <- polyfit(seq_along(y), y, k)
      y_pred <- polyval(fit, seq_along(y))
      y - y_pred
    }
    
    dat <- app(dat, detrend_pracma, k = k)
    # dat <- c(dat[[1]],dat)
  }
  
  
  if(isTRUE(anom) & length(ref_period == 2)){
    
    ref_mean <- dat %>% 
      subset(which(year(time_dat) >= ref_period[1] & year(time_dat) <= ref_period[2])) %>%
      app("mean", na.rm =T)
    
    if(is.function(conversion_fun)){
      ref_mean <- conversion_fun(ref_mean)    
    }
  
    
  }else if(isFALSE(anom)){
    
  } else if(isTRUE(anom) & is.null(ref_period)){
    stop("Please, a reference period is required. eg: ref_period = c(1981,2010)")
    
  }else{
    stop("Please, the reference period is as follows: ref_period = c(1981,2010)")
  }

  
  ts_nc <-tibble(id = seq_along(time(dat)), # id to subset
                 time = time_dat) 
  
  # Computing sd and mean for a counterfactual world -------------------------------
  
  # periods available
  yr_split <- analogs_subperiods$period %>% 
    unique() 
  
  sim_bs_l <- list()
  for(ii in seq_along(yr_split)){
    # selecting counterfactual analogs
    dat_subset <- inner_join(analogs_subperiods, ts_nc, by = "time") %>%
    filter(str_detect(period, yr_split[ii]))
    
      bootstrap <- function(x){
      cases <- dat_subset %>% group_by(time_obj) %>%
        slice_sample(n = 1, replace = replace)
    
      # Computing mean
      sim_event <- dat[[cases$id]] %>% # selecting VOI analogs
        app(event_fun) # conversion fun
  
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
