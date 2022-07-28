library(ecmwfr)
library(terra)
library(tidyverse)
library(lubridate)


download_era5 <- function(dataset = "ERA5_Land",
                          variable = "2m_temperature",
                          level = NULL,
                          event_dates = event_dates_,
                          time_window = 7,
                          year_range = c(1950,2022),
                          grid = 1,
                          area = c(-10,5,35,45),
                          time = time_,
                          out_dir = "data/",
                          hour2day_fun = "mean",
                          rm_hourly_files = T){
  
  area <- c(area[4],area[1],area[3],area[2])
  area <- paste(area,collapse = "/")
  timesteps <- paste(str_pad(time,width = 2,side = "left",pad = 0),collapse = "/")
  grid <- paste(grid,grid,sep = "/")
  
  if(is.null(level)){
    type = "single"
  }else{
    type = "pressure"
  }
  
  
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
  
  time_window_an <-  sapply(yr_seq, FUN = .time_windows) %>% 
    mdy() %>% 
    as_tibble() %>% 
    setNames("Dates") %>%
    separate(col = "Dates", into = c("yr","mo","dy"), sep = "-") %>%
    select(-yr) %>%
    distinct(mo,dy)
  
  
  mo2download <- pull(time_window_an,mo) %>% unique()
  
  for (oo in seq_along(mo2download)) {
  
  dy2download <- filter(time_window_an, mo == mo2download[oo] ) %>%
    pull(dy)
  
    for (yy in yr_seq) {
     if(dataset == "ERA5"){
      if(yy < 1959){
        product_type <- paste0('reanalysis-era5-',type,'-levels-preliminary-back-extension')
      } else {
        product_type <- paste0('reanalysis-era5-',type,'-levels')
      }
      
      if(type == "pressure"){
        request <- list(
          product_type = "reanalysis",
          format = "netcdf",
          variable = variable,
          grid = grid,
          pressure_level = as.character(level),
          year = as.character(yy),
          month = str_pad(mo2download, width = 2,side = "left","0"),
          day = str_pad(dy2download, width = 2,side = "left","0"),
          time = timesteps,
          area = area,
          dataset_short_name = product_type,
          target = paste0(variable,level,"_",yy,
                          str_pad(mo2download[oo], width = 2,side = "left","0"),
                          "_",str_replace_all(timesteps,pattern = "/",replacement = ""),
                          "_", dataset,".nc")
        )
      } else if(type == "single"){
        request <- list(
          product_type = "reanalysis",
          format = "netcdf",
          variable = variable,
          grid = grid,
          year = as.character(yy),
          month = str_pad(mo2download, width = 2,side = "left","0"),
          day = str_pad(dy2download, width = 2,side = "left","0"),
          time = timesteps,
          area = area,
          dataset_short_name = product_type,
          target = paste0(variable,"_",yy,
                          str_pad(mo2download[oo], width = 2,side = "left","0"),
                          "_",str_replace_all(timesteps,pattern = "/",replacement = ""),
                          "_", dataset,".nc")
        )
      }
     }else if(dataset == "ERA5_Land"){
       product_type <- 'reanalysis-era5-land'
       request <- list(
         format = "netcdf",
         variable = variable,
         grid = grid,
         year = as.character(yy),
         month = str_pad(mo2download, width = 2,side = "left","0"),
         day = str_pad(dy2download, width = 2,side = "left","0"),
         time = timesteps,
         area = area,
         dataset_short_name = product_type,
         target = paste0(variable,"_",yy,
                         str_pad(mo2download[oo], width = 2,side = "left","0"),
                         "_",str_replace_all(timesteps,pattern = "/",replacement = ""),
                         "_", dataset,".nc")
       )
     } else{
       stop("Choose only between ERA5 or ERA5_Land")
     }
      
    file <- wf_request(user     = "8933",   # user ID (for authentification)
                       request  = request,  # the request
                       transfer = TRUE,     # download the file
                       path     = out_dir)     # store data in current working directory
   
    if(!is.null(hour2day_fun)){
      if(hour2day_fun == "mean"){
        
        x <- rast(file)
        varname <- varnames(x)
        time <- as_date(time(x))
        x <- x %>% tapp(as.factor(time),"mean")
        time(x) <- unique(time)
        
        } else if(hour2day_fun == "max"){
          
          x <- rast(file)
          varname <- varnames(x)
          time <- as_date(time(x))
          x <- x %>% tapp(as.factor(time),"max")
          time(x) <- unique(time)
          
        } else if(hour2day_fun == "min"){
          
          x <- rast(file)
          varname <- varnames(x)
          time <- as_date(time(x))
          x <- x %>% tapp(as.factor(time),"min")
          time(x) <- unique(time)
        } else if(hour2day_fun == "sum"){
          
          x <- rast(file)
          varname <- varnames(x)
          time <- as_date(time(x))
          x <- x %>% tapp(as.factor(time),"sum")
          time(x) <- unique(time)
        } else{
          stop("For now no more options available than mean, max, min and sum")
        }
        
      
    }
    if(is.null(level)){
      writeCDF(x,
               filename = paste0(out_dir,variable,"_",yy,
                                 str_pad(unique(mo2download), width = 2,side = "left","0"),
                                 "_",str_replace_all(timesteps,pattern = "/",replacement = ""),
                                 "_", dataset,"_daily.nc"),
               varname = varname, 
               zname = "time", 
               prec = "float", 
               overwrite = T,
               compression = 3)
    } else{
      writeCDF(x,
               filename = paste0(out_dir,variable,level,"_",yy,
                                 str_pad(unique(mo2download), width = 2,side = "left","0"),
                                 "_",str_replace_all(timesteps,pattern = "/",replacement = ""),
                                 "_", dataset,"_daily.nc"),
               varname = varname, zname = "time", 
               prec = "float", 
               overwrite = T,compression = 3)
    }
    
    if(rm_hourly_files == T){
      if (file.exists(file)) {
        #Delete file if it exists
        file.remove(file)
      }
  }
  
  # if(concatenate == T){
  #   list.files(path = out_dir,pattern = )
  # }
  }
}

}
