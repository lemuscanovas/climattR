library(tidyverse)
library(lubridate)
library(furrr)
library(terra)
library(scales)

bs_analogs<- function(x, analogs_subperiods,
                      detrend = F,k=2, n = 1000, replace = T, anom = NULL, 
                      ref_period = NULL, hour2day_fun = "mean",
                      event_fun,conversion_fun,cl = 1 ){
  
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
      hour2day_FUN <- match.fun(FUN = hour2day_fun, na.rm = T)
      dat <- dat %>%
        group_by(time) %>%
        summarise(var = hour2day_FUN(var)) %>%
        ungroup()
    }
    dat <- dat %>%
      group_by(time) %>%
      summarise(var = mean(var, na.rm = T)) %>%
      ungroup()
  }    
  
  
  if(isTRUE(detrend)){
    detrend_pracma <- function(y, k) {
      fit <- polyfit(seq_along(y), y, k)
      y_pred <- polyval(fit, seq_along(y))
      y - y_pred
    }
    dat <- dat %>% na.omit() %>% mutate(var = detrend_pracma(var, k = k))
  }
  
  if(isTRUE(anom) & length(ref_period == 2)){
  ref_mean <- dat %>% 
    filter(year(time) >= ref_period[1], year(time) <= ref_period[2]) %>%
    summarise(ref_mean = mean(var, na.rm = T),
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
    unique() 
  
  # selecting n analogs per period 
  sim_bs_l <- list()
  for(ii in seq_along(yr_split)){
    
  dat_subset <- inner_join(analogs_subperiods, dat, by = "time") %>%
    filter(str_detect(period, (yr_split %>% str_sub(1,4))[ii])) %>% 
    relocate(period, .after = time)
  
  bootstrap <- function(x){
    cases <- dat_subset %>% group_by(time_obj) %>%
      slice_sample(n = 1,replace = replace) %>%
      mutate(sim = paste0("sim",x))
  }
  
  message(paste0("bootstrapping the event ", n," times in the ", yr_split[ii], ":"))
  
  sim_bs_l[[ii]] <- pbapply::pblapply(1:n, bootstrap, cl = cl) %>% 
                bind_rows()

  }
  
  sim_bs <- bind_rows(sim_bs_l)
 
  
  # COMPUTING STATISTICAL SIGN OF BOOTSTRAP  --------------------------------
  
  .range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  
  
  join_bootstraps <- sim_bs %>% 
    dplyr::select(sim,time_obj, dist, var,period) %>%
    # group_by(time_obj) %>%
    # mutate(dist = .range01(dist,na.rm = T)) %>%
    group_by(time_obj,sim,period) %>%
    summarise(var = ifelse(exists("conversion_FUN"), conversion_FUN(var),var),
           var = ifelse(exists("ref_mean"), var - pull(ref_mean,1),var),
           dist = mean(dist),
           var = var,.groups = "drop") %>%
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
              .groups = "drop") %>%
    ungroup()
  
  # mitjana anomalia de cada dia analeg
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
    summarise(value = ifelse(exists("conversion_FUN"), conversion_FUN(var),var),
           value = ifelse(exists("ref_mean"), value - pull(ref_mean,1),value),.groups = "drop") %>%
    ungroup() %>%
    rename(time_obj = time) 
  

  return(list(observed = dat_days_event,
              bootstrap_simulation = join_bootstraps,
              summary_bs = summary1_bootstrap_daily))
         
}
