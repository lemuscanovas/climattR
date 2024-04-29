
# Downscaling spatial analogs ---------------------------------------------


library(terra)
library(elevatr)
library(rnaturalearth)
library(gstat)   # The most popular R-Package for Kriging (imho)
library(automap) # Automatize some (or all) parts of the gstat-workflow 
library(sf)
library(gam)

# Downscaling validation event --------------------------------------------

spat_down <- function(reconstruction, disagg){

  # Reading and preparing data ------------------------------------------------------------
  # Reading analogs dates and VOI nc ----------------------------------------
  
  sample <- reconstruction[[1]][[1]] %>%
    setNames("z")
  # plot(sample)
  
  dem <- get_elev_raster(locations = sample,
                              prj = crs(sample),
                              z = 7) %>% 
    rast() 
  dem[dem < 0] <- NA
  
  x_coarse <- dem %>%
    project(sample)
  
  x_fine <- dem %>%
    project(disagg(sample,disagg))
  # 
  # plot(x_coarse)
  # plot(x_fine)
  
  ## Downscaling reconstructed periods
  periods <- names(reconstruction)

  .downscaling_dates <- function(ii){
    nn <- periods[ii]
    dd <- reconstruction[[ii]] %>% app("mean")

    day_sf <- as.points(dd) %>% st_as_sf() %>%
      cbind(st_coordinates(.))
    
    day_sf_ <- bind_cols(day_sf,
                        terra::extract(x_coarse,
                                       vect(day_sf)) %>% select(2) %>% 
                          setNames("elev")) %>%
      filter(!is.na(elev)) %>%
      rename("z"= 1)
    
    target <- as.points(x_fine) %>% st_as_sf() %>%
      cbind(st_coordinates(.)) %>% rename("elev" = 1)
    
    
    regression <- gam(z~s(X)+s(Y)+s(elev),
                      data = as.data.frame(day_sf_) %>% select(-ncol(.)))
    
    step_regression <- step.Gam(regression,trace = F,
                                scope=list("X"=~1+X+s(X,4)+s(X,6)+s(X,12),
                                           "Y"=~1+Y+s(Y,4)+s(Y,6)+s(Y,12),
                                           "elev"=~1+elev+s(elev,4)+s(elev,6)+s(elev,12))
    )
    
    xx <- formula(step_regression)
    
    v_mod_OK <- automap::autofitVariogram(xx, as(day_sf_, "Spatial"))$var_model
    
    mod <- krige(
      formula=formula(step_regression),
      locations=as_Spatial(day_sf_),
      newdata=as_Spatial(target),
      model=v_mod_OK,
      debug.level = -1)
    
    mod_r <- mod %>% as.data.frame() %>% rast(type = "xyz")
    # plot(mod_r$var1.pred-273.15)
    # plot(mod_r$var1.var)
    down_var <- mod_r$var1.pred
    names(down_var) <- nn
    return(down_var)
  }
  
  
  res <- pblapply(1:nlyr(datevent),FUN = .downscaling_dates) %>% rast()
  return(res)
}

