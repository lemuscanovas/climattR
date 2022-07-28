library(tidyverse)
library(lubridate)
library(furrr)
library(terra)
library(scales)

bs_analogs<- function(x, analogs_subperiods,
                      detrend = F, n = 1000, replace = T, anom = NULL, 
                      ref_period = NULL, hour2day_fun = "mean",
                      event_fun,conversion_fun ){
  
  event_FUN <- match.fun(FUN = event_fun)

  if(is.function(conversion_fun)){
  conversion_FUN <- match.fun(FUN = conversion_fun) 
  }
  
  # Reading analogs dates and VOI nc ----------------------------------------
  dat <- x %>%
    rename("var" = 2)
    
  time_dat <- as_date(dat$time)
  
  
  if(is.Date(time_dat) == F){
    stop("time is not a date object. Please, check time values in your data")
  }
  
  # If hourly, converts to daily
  if(length(unique(dat$time)) != length(unique(time_dat))){
    if(hour2day_fun != "mean"){
      hour2day_FUN <- match.fun(FUN = hour2day_fun)
      dat <- dat %>%
        group_by(time) %>%
        summarise(var = hour2day_FUN(var)) %>%
        ungroup()
    }
    dat <- dat %>%
      group_by(time) %>%
      summarise(var = mean(var)) %>%
      ungroup()
  }    
  
  
  if(isTRUE(detrend)){
    # linear detrend
    dat_yr_mean <- dat %>% 
      group_by(year(time)) %>%
      summarise(var = mean(var))
    
    lyrs <- 1:nrow(dat_yr_mean)
    
    .remove_slope <- function(x) { 
      m <- lm(x ~ lyrs) 
      rm_slope <- residuals(m) + predict(m)[1] # elimina slope
      return(rm_slope)
    }
    
    # remove slope from yearly time series
    dat_yr_rm_slope <- lapply(dat_yr_mean, .remove_slope)$var
    # getting slope for each year
    slope_yr <- dat_yr_mean$var - dat_yr_rm_slope
    # remove slope from daily ts
    yrs <- unique(year(time_dat))
    ind <- seq_along(yrs)
    
    .detrend_fun <-  function(ii){
      ind_yr <- which(year(time_dat)==yrs[ii])
      selyr_dat <- slice(dat, ind_yr)
      selyr_dat_slope <- slope_yr[ii]
      selyr_dat_detrended <- mutate(selyr_dat, var = var - selyr_dat_slope)
      
      return(selyr_dat_detrended)
    }
    
    # detrended
    dat <- lapply(ind, .detrend_fun) %>% bind_rows()
  }
  
  if(isTRUE(anom) & length(ref_period == 2)){
  ref_mean <- dat %>% 
    filter(year(time) >= ref_period[1], year(time) <= ref_period[2]) %>%
    summarise(ref_mean = mean(var),
              ref_mean = ifelse(is.function(conversion_fun), conversion_FUN(ref_mean),ref_mean))
  

    }else if(isFALSE(anom)){
      
    } else if(isTRUE(anom) & is.null(ref_period)){
    stop("Please, a reference period is required. eg: ref_period = c(1981,2010)")
    
  }else{
    stop("Please, the reference period is as follows: ref_period = c(1981,2010)")
      }

  # Computing sd and mean for a counterfactual world -------------------------------
  
  # to select first period (factual)
  yr_split <- analogs_subperiods$period %>% 
    str_sub(start = -4,end = -1) %>% 
    as.numeric() %>% 
    min %>% 
    as.character()
  
  # selecting counterfactual analogs (past)
  counter <- inner_join(analogs_subperiods, dat, by = "time") %>%
    filter(str_detect(period, yr_split)) %>% 
    relocate(period, .after = time)
  
  
  bootstrap_counter <- function(x){
    cases <- counter %>% group_by(time_obj) %>%
      slice_sample(n = 1,replace = replace) %>%
      mutate(sim = paste0("sim",x))
  }
  
  message(paste0("bootstrapping the event ", n," times in a counterfactual world:"))
  
  sim_bs_counter <- pbapply::pblapply(1:n, bootstrap_counter) %>% 
                bind_rows()

  
  # FACTUAL CALC (PRESENT)-----------------------------------------------------
  
  yr_split <- analogs_subperiods$period %>% 
    str_sub(start = -4,end = -1) %>% 
    as.numeric() %>% 
    min %>% sum(1) %>%
    as.character()
  
  # selecting counterfactual analogs
  factual <- inner_join(analogs_subperiods, dat, by = "time") %>%
    filter(str_detect(period, yr_split)) %>% 
    relocate(period, .after = time)
  
  
  bootstrap_factual <- function(x){
    cases <- factual %>% group_by(time_obj) %>%
      slice_sample(n = 1,replace = replace) %>%
      mutate(sim = paste0("sim",x))
  }
  
  message(paste0("bootstrapping the event ", n," times in a factual world:"))
  
  sim_bs_factual <- pbapply::pblapply(1:n, bootstrap_factual) %>% 
                bind_rows()
  
  
  
  
  # COMPUTING STATISTICAL SIGN OF BOOTSTRAP  --------------------------------
  
  .range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  
  
  join_bootstraps <- sim_bs_counter %>% 
    dplyr::select(sim,time_obj, dist, var,period) %>%
   bind_rows(
      dplyr::select(sim_bs_factual, 
               sim,time_obj, dist, var, period)) %>%
    group_by(time_obj) %>%
    mutate(dist = .range01(dist,na.rm = T)) %>%
    group_by(time_obj,sim,period) %>%
    mutate(var = ifelse(exists("conversion_FUN"), conversion_FUN(var),var),
           var = ifelse(exists("ref_mean"), var - pull(ref_mean,1),var)) %>%
    ungroup()
  
  # mean and sd of RMSD or Euclid dist, 
  # mean, median and sd of var
  summary1_bootstrap_daily <- join_bootstraps %>%
    group_by(time_obj,period) %>%
    summarise(dist_mean = mean(dist),
              dist_sd = sd(dist),
              var_mean = mean(var),
              var_median = median(var),
              var_sd = sd(var),
              .groups = "drop")
  
  
  # ttest of metric distance and var
  # summary2_bootstrap_daily <- join_bootstraps %>%
  #   pivot_wider(
  #     names_from = period,
  #     values_from = c(dist, var),values_fn = list) %>%
  #   unnest(cols = 2:5) %>%
  #   group_by(time_obj) %>%
  #   summarise(ttest_dist = t.test(.[[2]],.[[3]])$p.value,
  #             ttest_mean = t.test(.[[4]],.[[5]])$p.value)
  # 
  # summary2_bootstrap_event <- join_bootstraps %>%
  #   pivot_wider(
  #     names_from = period,
  #     values_from = c(dist, var),values_fn = list) %>%
  #   unnest(cols = 2:5) %>%
  #   summarise(ttest_dist = t.test(.[[2]],.[[3]])$p.value,
  #             ttest_mean = t.test(.[[4]],.[[5]])$p.value)
  

  dat_days_event <- dat %>% 
    filter(time %in% unique(analogs_subperiods$time_obj)) %>%
    group_by(time) %>%
    mutate(value = ifelse(exists("conversion_FUN"), conversion_FUN(var),var),
           value = ifelse(exists("ref_mean"), value - pull(ref_mean,1),value)) %>%
    ungroup() %>%
    rename(time_obj = time) %>%
    select(-var)
    
  return(list(observed = dat_days_event,
              bootstrap_simulation = join_bootstraps,
              summary_bs = summary1_bootstrap_daily))
         
}
