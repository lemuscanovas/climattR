
# Downscaling spatial analogs ---------------------------------------------


library(terra)
library(elevatr)
library(rnaturalearth)
library(gstat)   # The most popular R-Package for Kriging (imho)
library(automap) # Automatize some (or all) parts of the gstat-workflow 
library(sf)
library(gam)
sf_ob <- ne_countries(country = c("portugal","spain","france"),returnclass = "sf")
extent <- c(-3,4,41,44)
disagg = 2

t2m2 <- t2m %>% crop(extent)
# Downscaling validation event --------------------------------------------

spat_down <- function(x, event_dates, disagg){

  # Reading and preparing data ------------------------------------------------------------
  event_dates= seq(as.Date("2022-06-12"),as.Date("2022-06-21"), by = "day")
  
  t2m <- list.files(path = "inst/testdata/",pattern = "2m_t",full.names = T) %>% rast()
  datevent <- tidync_attr(x = t2m,level = NULL,detrend = F,scale = F,
                        extent = NULL, aggregate = NULL,rotate = F,
                        event_dates =  event_dates,time_window = 0,save = F )$event
  
  
  sample <- datevent[[1]] %>%
    setNames("z")
  # plot(sample)
  
  dem <- get_elev_raster(locations = raster::raster(sample),
                              prj = crs(sample),
                              z = 7) %>% 
    rast() 
  
  x_coarse <- dem %>%
    project(sample)
  
  x_fine <- dem %>%
    project(disagg(sample,disagg))
  # 
  # plot(x_coarse)
  # plot(x_fine)
  
  ## Downscaling event
  dates <- as_date(time(datevent))
  noms <- names(datevent)
  
  res <- pblapply(1:nlyr(datevent),FUN = .downscaling_dates) %>% rast()
  return(res)
}

# Downscaling -------------------------------------------------------------
.downscaling_dates <- function(ii){
  dd <- dates[ii]
  nn <- noms[ii] 
  x <- datevent[[ii]]
  names(x) <- "z"
  day_sf <- as.points(x) %>% st_as_sf() %>%
    cbind(st_coordinates(.))
  
  day_sf <- bind_cols(day_sf,
                      terra::extract(x_coarse,
                                     vect(day_sf)) %>% select(2) %>% 
                        setNames("elev")) %>%
    filter(!is.na(elev))
  
  target <- as.points(x_fine) %>% st_as_sf() %>%
    cbind(st_coordinates(.)) %>% rename("elev" = 1)
  
  
  regression <- gam(z~s(X)+s(Y)+s(elev),
                    data = as.data.frame(day_sf))
  
  step_regression <- step.Gam(regression,trace = F,
                              scope=list("X"=~1+X+s(X,4)+s(X,6)+s(X,12),
                                         "Y"=~1+Y+s(Y,4)+s(Y,6)+s(Y,12),
                                         "elev"=~1+elev+s(elev,4)+s(elev,6)+s(elev,12))
  )
  
  xx <- formula(step_regression)
  
  v_mod_OK <- automap::autofitVariogram(xx, as(day_sf, "Spatial"))$var_model
  
  mod <- krige(
    formula=formula(step_regression),
    locations=as_Spatial(day_sf),
    newdata=as_Spatial(target),
    model=v_mod_OK,
    debug.level = 0)
  
  mod_r <- mod %>% as.data.frame() %>% rast(type = "xyz")
  # plot(mod_r$var1.pred-273.15)
  # plot(mod_r$var1.var)
  down_var <- mod_r$var1.pred
  time(down_var) <- dd
  names(down_var) <- nn
  return(down_var)
}
