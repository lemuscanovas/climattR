#' Downscale Temperature Using Spatial Analogs
#'
#' This function downscales temperature data using spatial analogs. It integrates various
#' geographic and statistical methodologies to refine coarse-resolution temperature data
#' to a finer resolution based on digital elevation models (DEM) and other geographical
#' transformations. The function uses Kriging for spatial interpolation.
#'
#' @param x A list of raster layers representing different periods or scenarios.
#' @param dem An optional digital elevation model raster. If not provided, it will be automatically
#'        downloaded using the `elevatr` package.
#' @param disagg A numeric value indicating the disaggregation factor for downscaling.
#' @param z An integer specifying the zoom level for DEM retrieval if `dem` is not provided.
#'
#' @return A raster object representing the downscaled temperature data across different periods.
#'
#' @import terra
#' @import elevatr
#' @import rnaturalearth
#' @import gstat
#' @import automap
#' @import sf
#' @import gam
#' @examples
#' # Assuming 'temp_rasters' is a list of raster layers:
#' # temp_rasters <- list(rast1, rast2, rast3)
#' # dem_raster <- rast("path_to_dem.tif")
#' # downscaled_temps <- downscale_temperature(temp_rasters, dem = dem_raster, disagg = 2, z = 6)
#' @export
downscale_temperature <- function(x, dem = NULL, disagg, z = 6) {

  # Reading and preparing data ------------------------------------------------------------
  sample <- x[[1]][[1]] %>%
    setNames("z")

  if(is.null(dem)){
  dem <- elevatr::get_elev_raster(locations = sample,
                              prj = crs(sample),
                              z = z) %>% 
    rast() 
  }
  
  dem[dem < 0] <- NA
  
  x_coarse <- dem %>%
    project(sample)
  
  x_fine <- dem %>%
    project(disagg(sample,disagg))

  
  ## Downscaling reconstructed periods
  periods <- names(x)

  .downscaling_reconstructions <- function(ii){
    nn <- periods[ii]
    dd <- x[[ii]] %>% app("mean")

    day_sf <- as.points(dd) %>% sf::st_as_sf() %>%
      cbind(sf::st_coordinates(.))
    
    day_sf_ <- bind_cols(day_sf,
                        terra::extract(x_coarse,
                                       vect(day_sf)) %>% dplyr::select(2) %>% 
                          setNames("elev")) %>%
      filter(!is.na(elev)) %>%
      rename("z"= 1)
    
    target <- as.points(x_fine) %>% sf::st_as_sf() %>%
      cbind(sf::st_coordinates(.)) %>% rename("elev" = 1)
    
    
    regression <- gam::gam(z~s(X)+s(Y)+s(elev),
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

    down_var <- mod_r$var1.pred
    names(down_var) <- nn
    message(paste0("Spatial downscaling provided for reconstructed period: ", nn))
    
    return(down_var)
  }
  
  
  res <- lapply(seq_along(periods),FUN = .downscaling_reconstructions) %>% rast()
  return(res)
}

