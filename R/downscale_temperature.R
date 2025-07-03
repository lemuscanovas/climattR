#' Downscale Temperature Using Spatial Analogs
#'
#' This function downscales temperature data using spatial analogs. It integrates various
#' geographic and statistical methodologies to refine coarse-resolution temperature data
#' to a finer resolution based on digital elevation models (DEM) and other geographical
#' transformations. The function uses Kriging for spatial interpolation.
#'
#' @param x A list of raster layers or a SpatRaster object representing different periods or scenarios.
#' @param analogs A list specifying the time periods for the analogs used in downscaling.
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
#' @import giscoR
#' @examples
#' \dontrun{
#' # Assuming 'temp_rasters' is a list of raster layers:
#' # temp_rasters <- list(rast1, rast2, rast3)
#' # analogs <- list(time = c("2001-01-01", "2001-01-02", "2001-01-03"))
#' # dem_raster <- rast("path_to_dem.tif")
#' # downscaled_temps <- downscale_temperature(temp_rasters, analogs = analogs, dem = dem_raster, disagg = 2, z = 6)}
#' @export

downscale_temperature <- function(x, analogs= analogs, dem = NULL, disagg, z = 6) {

  
  # Initialize raster
if (inherits(x, "SpatRaster")) {
    dat <- x
  } else {
    dat <- rast(x)
  }
  
  proj_ <- crs(dat)
  # slect analog fields
  dat <- dat[[which(terra::time(dat) %in% unique(analogs$time))]]
  world <- giscoR::gisco_get_coastallines(resolution = 10)
  
  if(is.null(dem)){
  dem <- elevatr::get_elev_raster(locations = dat[[1]],
                              prj = crs(dat[[1]]),
                              z = z) %>% 
    rast() 
  }
  
  dem <- dem %>% mask(world)
  
  x_coarse <- dem %>%
    project(dat[[1]])
  
  x_fine <- dem %>%
    project(disagg(dat[[1]],disagg))

  
  ## Downscaling reconstructed periods
  times <- time(dat)

  .downscaling_reconstructions <- function(ii){
    nn <- times[ii]
    dd <- dat[[ii]]

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
    
    
    regression <- lm(z~X+Y+elev,
                      data = as.data.frame(day_sf_) %>% select(-ncol(.)))
    
    step_regression <- step(regression,direction = "backward",trace = 0)
    
    xx <- formula(step_regression)
    
    v_mod_OK <- automap::autofitVariogram(xx, as(day_sf_, "Spatial"))$var_model
    
    mod <- krige(
      formula=formula(step_regression),
      locations=sf::as_Spatial(day_sf_),
      newdata=sf::as_Spatial(target),
      model=v_mod_OK,
      debug.level = 0)
    
    mod_r <- mod %>% as.data.frame() %>% rast(type = "xyz")

    down_var <- mod_r$var1.pred
    time(down_var) <- nn
    # message(paste0("Spatial downscaling provided for reconstructed period: ", nn))
    
    return(down_var)
  }
  
  
  res <- pblapply(seq_along(times),FUN = .downscaling_reconstructions) %>% rast()
  crs(res) <-proj_
  return(res)
}

