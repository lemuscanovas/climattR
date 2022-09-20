library(metR)
library(tidyverse)
library(lubridate)
library(terra)

analogs_searcher <- function(ts_wo_event, event, n = 20, split_year, metric = "rmsd"){
  
  if(metric == "rmsd"){
    FUN <- function(a, b) sqrt(sum((a - b)^2)/length(a))
  } else if(metric  == "euclidean"){
    FUN <- function(a, b) sqrt(sum((a - b)^2))
  }else{
    stop("Choose between: rmsd or euclidean")
  }
  
  ts_wo_event_df <-  ts_wo_event %>% as.data.frame(xy = T) %>% as_tibble() %>%
    pivot_longer(names_to = "time",values_to = "var_hist", 3:ncol(.)) %>%
    mutate(time = ymd(str_replace(time,"X","")),
           period = ifelse(year(time) > 1985,
                           str_c(split_year+1,last(year(time)),sep = "-"),
                           str_c(first(year(time)),split_year,sep ="-")))
  
  event_df <- event %>% as.data.frame(xy = T) %>% as_tibble() %>%
    pivot_longer(names_to = "time",values_to = "var_obj", 3:ncol(.)) %>%
    mutate(time = ymd(str_replace(time,"X","")))
  
  
  time_event = unique(event_df$time)
  
  # FUN <- match.fun(metric)
  .range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  
  .analog_extraction_all_period <- function(s){
    # Dates anàlogues a z500
    event_df_ii <- filter(event_df, time == time_event[s]) %>%
      select(-time)
    
    analog_dates <- inner_join(ts_wo_event_df,
                               event_df_ii, 
                               by = c("x","y")) %>%
      group_by(time) %>%
      summarise(dist = .range01(FUN(var_hist,
                           var_obj)),
                .groups = "drop") %>%
      ungroup() %>%
      arrange(dist) %>%
      slice(1:n) %>%
      mutate(time_obj = as_date(time_event[s])) %>%
      select(time_obj, time, dist)  
  }
  
  .analog_extraction_subperiods <- function(s){
    # Dates anàlogues a z500
    event_df_ii <- filter(event_df, time == time_event[s]) %>%
      select(-time)
    
    analog_dates <- inner_join(ts_wo_event_df,
                               event_df_ii, 
                               by = c("x","y")) %>%
      group_by(time,period) %>%
      summarise(dist = .range01(FUN(var_hist,
                           var_obj)),
                .groups = "drop") %>%
      group_by(period) %>%
      arrange(dist) %>%
      slice(1:n) %>%
      ungroup() %>%
      mutate(time_obj = as_date(time_event[s])) %>%
      select(time_obj, time, dist, period)  
  }
  
  analogs_full_period <- lapply(seq_along(time_event), .analog_extraction_all_period) %>% bind_rows()
  analogs_subperiods <- lapply(seq_along(time_event), .analog_extraction_subperiods) %>% bind_rows()
  
  return(list(analogs_full_period = analogs_full_period,
              analogs_subperiods = analogs_subperiods))
}

