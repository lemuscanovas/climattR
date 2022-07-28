library(tidyverse)
library(lubridate)

bs_spatanalogs<- function(x, analogs_subperiods,detrend = F,
                                 n = 1000, replace = T,
                                 event_fun,
                                hour2day_fun = "mean",
                                conversion_fun = NULL,
                                 save = F){
  
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
  
  ts_nc <-tibble(id = seq_along(time(dat)), # id to subset
                 time = time_dat) 
  
  # Computing sd and mean for a counterfactual world -------------------------------
  
  # to select first period (counterfactual)
  yr_split <- analogs_subperiods$period %>% 
    str_sub(start = -4,end = -1) %>% 
    as.numeric() %>% 
    min %>% 
    as.character()
  
  # selecting counterfactual analogs
  counter <- inner_join(analogs_subperiods, ts_nc, by = "time") %>%
    filter(str_detect(period, yr_split))
  
  
  bootstrap_counter <- function(x){
    cases <- counter %>% group_by(time_obj) %>%
      slice_sample(n = 1, replace = replace)
  
    # Computing mean
    sim_event <- dat[[cases$id]] %>% # selecting VOI analogs
      app(event_fun) # conversion fun

  }

message(paste0("bootstrapping the event ", n," times in a counterfactual world:"))
  
system.time(sim_bs_counter <- pbapply::pblapply(1:n, bootstrap_counter) %>% 
              rast())
names(sim_bs_counter) <- paste0("sim",1:n)


# FACTUAL CALC -----------------------------------------------------
  
  yr_split <- analogs_subperiods$period %>% 
    str_sub(start = -4,end = -1) %>% 
    as.numeric() %>% 
    min %>% sum(1) %>%
    as.character()
  
  # selecting factual analogs
  factual <- inner_join(analogs_subperiods, ts_nc, by = "time") %>%
    filter(str_detect(period, yr_split))
  
  bootstrap_factual <- function(x){
    cases <- factual %>% group_by(time_obj) %>%
      slice_sample(n = 1, replace = replace)
    
    # Computing mean
    sim_event <- dat[[cases$id]] %>% # selecting VOI analogs
      app(event_fun) # conversion fun
    
  }

message(paste0("bootstrapping the event ", n," times in a factual world:"))
system.time(sim_bs_factual <- pbapply::pblapply(1:n, bootstrap_factual) %>% 
              rast())
names(sim_bs_factual) <- paste0("sim",1:n)




# COMPUTING STATISTICAL SIGN OF BOOTSTRAP  --------------------------------

ttest_bootstrap <- as.data.frame(sim_bs_factual, xy = T) %>% 
  tidyr::as_tibble() %>% 
  pivot_longer(names_to = "sim_f",values_to = "value_f",3:ncol(.)) %>%
  bind_cols(
    select(as.data.frame(sim_bs_counter, xy = T) %>% 
           tidyr::as_tibble() %>% 
           pivot_longer(names_to = "sim_cf",
                        values_to = "value_cf",
                        3:ncol(.)),
           value_cf)) %>%
  group_by(x,y) %>%
  summarise(mean_factual = mean(value_f),
            mean_counter= mean(value_cf),
            p.value = t.test(value_f,value_cf)$p.value,
            .groups = "drop") %>%
  rast(type="xyz") 

  if(!is.null(conversion_fun)){
  conversion_FUN <- match.fun(FUN = conversion_fun)
  a <- ttest_bootstrap[[1:2]] %>%
    app(conversion_FUN)
  }else{
    a <- ttest_bootstrap[[1:2]]
  }

b <- ttest_bootstrap[[1]] - ttest_bootstrap[[2]]
c <- ttest_bootstrap[[3]]

result_event_bootstrap <- c(a,b,c)
names(result_event_bootstrap)[3] <- "mean_dif"
  
  # Saving objects ----------------------------------------------------------
  if(isTRUE(save)){
    folder_name <- "02_bootstrap_spatanalogs"
    dir.create(file.path(folder_name), showWarnings = FALSE)      
    
    writeCDF(result_event_bootstrap,
             filename = paste0(folder_name,"/", varnames(dat),"_bootstrap_event.nc"),
             varname = varnames(dat),
             prec = "float", 
             overwrite = T,compression = 3)
    }
  
  return(result_event_bootstrap)
}