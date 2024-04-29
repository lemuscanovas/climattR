library(metR)
library(tidyverse)
library(lubridate)
library(terra)

analogs_searcher <- function(ts_wo_event, event, n = 20, 
                             periods = c(1951,1980,1991,2020), 
                             metric = "rmsd"){
  
  if(metric == "rmsd"){
    FUN <- function(a, b) sqrt(sum((a - b)^2)/length(a))
  } else if(metric  == "euclidean"){
    FUN <- function(a, b) sqrt(sum((a - b)^2))
  }else if(metric == "pearson"){
    FUN <- function(a, b) -(cor(a,b))
  }else{
    stop("Choose between: rmsd or euclidean")
  }
  
  ts_wo_event_df_l <- list()
  event_df_l <- list()
  
  for (var in 1:length(ts_wo_event)) {
    
    # event merging variables
    event_df_l[[var]] <- event[[var]] %>% 
      as.data.frame(xy = T) %>% 
      as_tibble() %>%
      pivot_longer(names_to = "time", values_to = "var_obj", 3:ncol(.)) %>%
      mutate(time = ymd(str_replace(time,"X","")),var_id = var) 
    
    # data outside event merging
    coords <- (event[[var]] %>% 
                 as.data.frame(xy = T))[,1:2]
    
    new_names <- str_c("X.",names(ts_wo_event[[var]]))
    
    names(ts_wo_event[[var]]) <- new_names
    
    ts_wo_event_tidy <- ts_wo_event[[var]] %>%
      as.data.frame(xy =T) %>% 
      as_tibble() %>% 
      filter(x %in% coords$x) %>% 
      select(-c(x,y))
    
    if(!is.null(periods)){
     
      split_year <- periods
      l = seq(1,length(split_year),2)
      
      splitted_data_l <- list()
      for (ff in l) {
        
        splitted_data_l[[ff]] <- bind_cols(coords, ts_wo_event_tidy) %>%
          pivot_longer(names_to = "time",values_to = "var_hist", 3:ncol(.)) %>%
          mutate(time = ymd(str_replace(time,"X",""))) %>% 
          filter(year(time) >= split_year[ff],year(time) <= split_year[ff+1]) %>%
          mutate(period = paste0(split_year[ff],"-",split_year[ff+1])) 
        
      }
      splitted_data <- bind_rows(splitted_data_l)
      
      ts_wo_event_df_l[[var]] <- splitted_data %>%
        mutate(var_id = var)
      
    }else{

    ts_wo_event_df_l[[var]] <- bind_cols(coords, ts_wo_event_tidy) %>%
      pivot_longer(names_to = "time",values_to = "var_hist", 3:ncol(.)) %>%
      mutate(time = ymd(str_replace(time,"X","")),var_id = var)
    }
  }
  
  ts_wo_event_df <- bind_rows(ts_wo_event_df_l)
  event_df <- bind_rows(event_df_l)
  
  time_event = unique(event_df$time)
  
  # FUN <- match.fun(metric)
  .range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  
  .analog_extraction_all_period <- function(s){
    # Dates anàlogues a z500
    event_df_ii <- filter(event_df, time == time_event[s]) %>%
      select(-time)
    
    analog_dates <- inner_join(ts_wo_event_df,
                               event_df_ii,
                               by = c("x","y", "var_id")) %>%
      group_by(time,var_id) %>%
      summarise(dist = FUN(var_hist, var_obj),.groups = "drop") %>%
      ungroup() %>%
      group_by(time) %>%
      summarise(dist = sum(dist), .groups = "drop") %>%
      ungroup() %>%
      arrange(dist) %>%
      slice(1:n) %>%
      mutate(time_obj = as_date(time_event[s])) %>%
      select(time_obj, time, dist)
  } 
  
  if(!is.null(split_year)){
  .analog_extraction_subperiods <- function(s){
    # Dates anàlogues a z500
    event_df_ii <- filter(event_df, time == time_event[s]) %>%
      select(-time)
    
    analog_dates <- inner_join(ts_wo_event_df, 
                               event_df_ii, 
                               by = c("x","y", "var_id")) %>%
      group_by(time,period, var_id) %>%
      summarise(dist = FUN(var_hist, var_obj), .groups = "drop") %>%
      ungroup() %>%
      group_by(period, time) %>%
      summarise(dist = sum(dist), .groups = "drop") %>%
      group_by(period) %>%
      arrange(dist) %>%
      slice(1:n) %>%
      ungroup() %>%
      mutate(time_obj = as_date(time_event[s])) %>%
      select(time_obj, time, dist, period)
    
  }
  analogs_subperiods <- lapply(seq_along(time_event), 
                               .analog_extraction_subperiods) %>% 
    bind_rows()
  analogs_full_period <- lapply(seq_along(time_event), 
                                .analog_extraction_all_period) %>%
    bind_rows()
  
  return(list(analogs_full_period = analogs_full_period,
              analogs_subperiods = analogs_subperiods))
  } else{
  analogs_full_period <- lapply(seq_along(time_event), 
                                .analog_extraction_all_period) %>% 
    bind_rows()
  
  return(list(analogs_full_period = analogs_full_period))
  }
}

