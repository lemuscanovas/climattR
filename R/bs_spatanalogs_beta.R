library(tidyverse)
library(lubridate)
library(pracma)
bs_spatanalogs<- function(x,sf=NULL,
                          analogs_subperiods,detrend = F,k = 2,
                                 n = 1000, replace = T,
                                 event_fun,
                                hour2day_fun = "mean",
                                conversion_fun = NULL,
                                anom = NULL, 
                                ref_period = NULL,
                                 save = F,
                          cl = 1){
  
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
    
    if(!is.null(sf)){
      ref_mean <- (dat %>% 
        subset(which(year(time_dat) >= ref_period[1] & year(time_dat) <= ref_period[2])) %>%
        app("mean") %>% terra::extract(vect(sf), fun = "mean", na.rm = T ))$mean
    }
      
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
  
system.time(sim_bs_counter <- pbapply::pblapply(1:n, bootstrap_counter, cl = cl) %>% 
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
system.time(sim_bs_factual <- pbapply::pblapply(1:n, bootstrap_factual,cl = cl) %>% 
              rast())
names(sim_bs_factual) <- paste0("sim",1:n)




# COMPUTING STATISTICAL SIGN OF BOOTSTRAP  --------------------------------

if(is.null(sf)){
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
    summarise(mean_factual = mean(value_f, na.rm = T),
              mean_counter= mean(value_cf, na.rm =T),
              p.value = t.test(value_f,value_cf)$p.value,
              .groups = "drop") %>%
    rast(type="xyz") 
  
    crs(ttest_bootstrap) <- crs(dat)
    
    if(!is.null(conversion_fun)){
    conversion_FUN <- match.fun(FUN = conversion_fun)
    a <- ttest_bootstrap[[1:2]] %>%
      app(conversion_FUN)
    }else{
      a <- ttest_bootstrap[[1:2]]
    }
  
  if(isTRUE(anom)){
    a <- a - project(ref_mean, a)
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
               filename = paste0(folder_name,"/", unique(varnames(dat)),"_bootstrap_event.nc"),
               varname = unique(varnames(dat)),
               prec = "float", 
               overwrite = T,compression = 3)
    }
  return(result_event_bootstrap)
}
 if(!is.null(sf)){
    
    if(class(sf)[1] == "sfc_POLYGON"){    
    counter_o <- terra::extract(sim_bs_counter,vect(sf), fun = mean, na.rm = T ) %>%
      t() %>% as.vector()
  
    factual_o  <- terra::extract(sim_bs_factual,vect(sf), fun = mean, na.rm = T ) %>%
      t() %>% as.vector()
    } else{
      counter_o <- terra::extract(sim_bs_counter,vect(sf), fun = NULL, na.rm = T ) %>%
        t() %>% as.vector()
      
      factual_o  <- terra::extract(sim_bs_factual,vect(sf), fun = NULL, na.rm = T ) %>%
        t() %>% as.vector()
    }
    
    
    periods <- analogs1$analogs_subperiods$period %>% unique
    
    join_bootstraps <- tibble(names(sim_bs_counter), counter_o[-1],factual_o[-1]) %>% 
      setNames(c("sim", periods[1], periods[2]))  %>%
      pivot_longer(2:3,names_to = "period",values_to = "var") %>%
      group_by(sim,period) %>%
      summarise(var = ifelse(exists("conversion_FUN"), conversion_FUN(var),var),
                var = ifelse(exists("ref_mean"), var - ref_mean,var),
                var = var,.groups = "drop") %>%
      ungroup()
    
    # mean and sd of RMSD or Euclid dist, 
    # mean, median and sd of var
    summary1_bootstrap_daily <- join_bootstraps %>%
      group_by(period) %>%
      summarise(#dist_mean = mean(dist),
                #dist_sd = sd(dist),
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
    
    
    dat_days_event <- dat[[which(as_date(time(x)) %in% unique(analogs_subperiods$time_obj))]] %>%
      terra::extract(vect(sf), fun = mean, na.rm = T ) %>%
      t() %>% as.data.frame() %>% as_tibble() %>% slice(-1) %>%
      mutate(time = unique(analogs_subperiods$time_obj)) %>%
      rename("var" = 1) %>% 
      group_by(time) %>%
      summarise(value = ifelse(exists("conversion_FUN"), conversion_FUN(var),var),
                value = ifelse(exists("ref_mean"), value - ref_mean,value),.groups = "drop") %>%
      ungroup() %>%
      rename(time_obj = time) 
    
    
    return(list(observed = dat_days_event,
                bootstrap_simulation = join_bootstraps,
                summary_bs = summary1_bootstrap_daily))
    
  }
}
