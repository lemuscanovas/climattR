
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
disagg = 10
# Downscaling validation event --------------------------------------------

spat_down <- function(x, event_dates, disagg)

# Reading and preparing data ------------------------------------------------------------
event_dates= seq(as.Date("2022-06-12"),as.Date("2022-06-21"), by = "day")

t2m <- list.files(path = "inst/testdata/",pattern = "2m_t",full.names = T) %>% rast()
datevent <- tidync_attr(x = t2m,level = NULL,detrend = F,scale = F,
                      extent = NULL, aggregate = NULL,rotate = F,
                      event_dates =  event_dates,time_window = 0,save = F )$event %>%
  mask(vect(sf_ob))


sample <- datevent[[1]] %>%
  setNames("z")
# plot(sample)

x_coarse <- get_elev_raster(locations = raster::raster(sample),
                            prj = crs(sample),
                            z = 7) %>% 
  rast() %>%
  project(sample) %>%
  mask(vect(ne_countries(continent = "europe",scale = 10, returnclass = "sf")))

x_fine <- get_elev_raster(locations = raster::raster(sample),
                          prj = crs(sample),
                          z = 7) %>% 
  rast() %>%
  project(disagg(sample,disagg)) %>%
  mask(vect(ne_countries(continent = "europe",scale = 10, returnclass = "sf")))
# 
# plot(x_coarse)
# plot(x_fine)

## Downscaling event
dates <- as_date(time(datevent))
noms <- names(datevent)


 # Downscaling -------------------------------------------------------------
downscaling_dates <- function(ii){
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
  down_var <- mod_r$var1.pred - 273.15
  time(down_var) <- dd
  names(down_var) <- nn
  return(down_var)
}

res <- pblapply(1:nlyr(datevent),FUN = downscaling_dates) %>% rast()

# Downscaling factual i counterfactual ------------------------------------

factual <- rast("01_counter_factual_spatanalogs/t2m_counterfactual_event_mean.nc") %>%
  crop(c(1,2,42.2,43)) %>%
  setNames("z")
plot(factual)

x_coarse <- get_elev_raster(locations = raster::raster(factual),
                prj = crs(factual),
                z = 7) %>% 
  rast() %>%
  project(factual) %>%
  mask(vect(ne_countries(continent = "europe",scale = 10, returnclass = "sf")))

x_fine <- get_elev_raster(locations = raster::raster(factual),
                          prj = crs(factual),
                          z = 7) %>% 
  rast() %>%
  project(disagg(factual,2)) %>%
  mask(vect(ne_countries(continent = "europe",scale = 10, returnclass = "sf")))

plot(x_coarse)
plot(x_fine)


# Downscaling -------------------------------------------------------------


library(gstat)   # The most popular R-Package for Kriging (imho)
library(automap) # Automatize some (or all) parts of the gstat-workflow 
library(sf)
library(gam)

factual_sf <- as.points(factual) %>% st_as_sf() %>%
  cbind(st_coordinates(.))
        
factual_sf <- bind_cols(factual_sf,
                        extract(x_coarse, 
                vect(factual_sf)) %>% select(2) %>% 
          setNames("elev"))

target <- as.points(x_fine) %>% st_as_sf() %>%
  cbind(st_coordinates(.)) %>% rename("elev" = 1)

regression <- lm(z~X+Y+elev, 
                   data = as.data.frame(factual_sf))
regression <- gam(z~s(X)+s(Y)+s(elev),
                 data = as.data.frame(factual_sf))

step_regression <- step(regression)
step_regression <- step.Gam(regression,
                            scope=list("X"=~1+X+s(X,4)+s(X,6)+s(X,12),
                                       "Y"=~1+Y+s(Y,4)+s(Y,6)+s(Y,12),
                                       "elev"=~1+elev+s(elev,4)+s(elev,6)+s(elev,12))
)

xx <- formula(step_regression)

v_mod_OK <- automap::autofitVariogram(xx, as(factual_sf, "Spatial"))$var_model
plot(automap::autofitVariogram(xx, as(factual_sf, "Spatial")))
n = round(nrow(factual_sf)/2)

mod <- krige(
  formula=formula(step_regression),
  locations=as_Spatial(factual_sf),
  newdata=as_Spatial(target),
  model=v_mod_OK)
mod_r <- mod %>% as.data.frame() %>% rast(type = "xyz")
plot(mod_r$var1.pred)

cv <- krige.cv(formula = formula(step_regression),
         locations =as_Spatial(factual_sf),
         v_mod_OK,
         nfold=n)

Metrics::rmse(cv$var1.pred,cv$observed)

mod_r <- mod %>% as.data.frame() %>% rast(type = "xyz") %>% 
  crop(c(-1,1,42,43))
plot(mod_r$var1.pred)
