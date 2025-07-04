# climattR <img src="man/figures/logo.png" align="right" width="140"/> 




# Rapid climate extreme event attribution for regional and local areas

[![CRAN status](https://www.r-pkg.org/badges/version/climattR)](https://cran.r-project.org/package=climattR) [![Project Status: Active -- The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) 
[![DOI](https://zenodo.org/badge/518908035.svg)](https://doi.org/10.5281/zenodo.15772289)

## Overview

**climattR** is a toolbox for performing rapid attribution of extreme weather events using the flow analogs approach. *This package is still under developement*

------------------------------------------------------------------------

Contents:

-   [Why this package](#why-this-package)
-   [How it works](#how-it-works)
    -   [Installation](#installation)
-   [Package citation](#package-citation)
-   [Contact](#contact)

------------------------------------------------------------------------

## Why this package {#why-this-package}

To date there is no R package available which have implemented the simple but effective flow analogues methodology to perform extreme weather event attribution. `climattR` contains a set of functions to 1) prepare data to employ flow analogues; 2) Find flow analogue days to a target date; 3) Reconstruct surface fields for extreme weather events in past and present (also future) climate conditions;

## How it works {#how-it-works}

### Installation {#installation}

``` r
# To install the latest version from Github:
# install.packages("remotes")
# remotes::install_github("lemuscanovas/climattR")

library(climattR)
```

### First steps. Rapid attribution of

Additional libraries

``` r
library(terra)
library(tidyterra)
library(tidyverse)
library(patchwork)
library(giscoR)

# terra::terraOptions(todisk = T)
borders <- gisco_get_coastallines()
```

``` r
# Loading daily Z500 and t2m data for 1950-2023 summers

Z500_file <-  system.file("extdata", "z500_0509_1950_2023_eu.nc", package = "climattR")
Tx_file <-  system.file("extdata", "tx_0608_1950_2023.nc", package = "climattR")

Z500 <- rast(Z500_file) 
Tx <- rast(Tx_file)
```

``` r
## Preparing data for analogs
dates <- seq(as_date("2023-07-12"), as_date("2023-07-18"),"day")

prepared_data <- prepare_data(x = Z500,
                       level = NULL, # the pressure level is already selected
                       event_dates =  dates,
                       time_window = 31)

extent <- ext(prepared_data$event)                
ggplot() +
  geom_spatraster(data = prepared_data$event - app(prepared_data$ts_wo_event, "mean"),
                  interpolate = T) +
  geom_spatraster_contour_text(data = prepared_data$event - app(prepared_data$ts_wo_event, "mean"),
                          breaks = seq(-200, 200, 40)) +
  geom_sf(data = borders, fill = "transparent", color = "black")+
  facet_wrap(~lyr, ncol = 3) +
  metR::scale_fill_divergent("m")+
  guides(fill = guide_colourbar(theme = theme(
         legend.key.width  = unit(8, "lines"),
         legend.key.height = unit(0.5, "lines"),
         legend.title.position = "top",
         legend.title = element_text(hjust = 0.5)))) + 
  scale_x_continuous(limits = c(extent[1],extent[2]), expand = c(0,0))+
  scale_y_continuous(limits = c(extent[3],extent[4]), expand = c(0,0))+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.7,0.15),
        legend.direction = "horizontal")
```
<img src="man/figures/event_z500.png" alt="" width="800"/>


``` r
## Searching analogs
analogs <- analogs_searcher(ts_wo_event = list(prepared_data$ts_wo_event),
                            event = list(prepared_data$event),
                            n = 20,
                            periods = c(1951,1980,1993,2022),
                            metric = "rmsd")

# Reporting the quality of the analogues       
ggplot(data = analogs$analogs_subperiods, aes(x = period, y = dist))+ 
geom_boxplot()+
labs(y = "RMSD", subtitle = "Analogs quality")+
theme_bw()+
theme(axis.title.x = element_blank())
```
<img src="man/figures/analogs_quality.png" alt="" width="300"/>

```r
# Reporting analogues frequency 

annual_counts <- analogs$analogs_full_period %>% mutate(yr = year(time)) %>% 
                 group_by(yr) %>% 
                 summarise(l = length(dist))
                 
ggplot(data = annual_counts, aes(x = yr, y = l))+ 
geom_col()+
labs(y = "Frequency", subtitle = "Analogs trend")+
theme_bw()+
theme(axis.title.x = element_blank())
```
<img src="man/figures/analogs_annual_freq.png" alt="" width="500"/>

```r
## Bootstrap heatwave Z500 reconstruction for a counterfactual/factual world
bs_sp_z500 <- bs_spatanalogs(x = Z500,
                                     analogs = analogs$analogs_subperiods,
                                     n = 1000,
                                     event_fun = "mean", 
                                     detrend = F,
                                     replace = T,
                                     anom = T,
                                     ref_period = c(1950,2022))
                                     
counterfactual <- ((bs_sp_z500$`1951-1980`) %>% app("mean"))
factual <- ((bs_sp_z500$`1993-2022`) %>% app("mean"))
dif <- factual - counterfactual
all <- c(counterfactual,factual,dif)
names(all) <- c(names(bs_sp_z500),"dif")
extent <- ext(all)                

a_z500 <- ggplot() +
  geom_spatraster(data = all[[-3]],interpolate = T) +
  geom_spatraster_contour_text(data = all[[-3]], breaks = seq(-200, 200, 20)) +
  geom_sf(data = borders, fill = "transparent", color = "black")+
  facet_wrap(~lyr, ncol = 3) +
  metR::scale_fill_divergent("m", na.value = "transparent")+
  guides(fill = guide_colourbar(theme = theme(
         legend.key.width  = unit(8, "lines"),
         legend.key.height = unit(0.5, "lines"),
         legend.title.position = "top",
         legend.title = element_text(hjust = 0.5)))) + 
  scale_x_continuous(limits = c(extent[1],extent[2]), expand = c(0,0))+
  scale_y_continuous(limits = c(extent[3],extent[4]), expand = c(0,0))+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")
        
b_z500 <- ggplot() +
  geom_spatraster(data = all[[3]],interpolate = T) +
  geom_spatraster_contour_text(data = all[[3]], breaks = seq(-20, 20, 5)) +
  geom_sf(data = borders, fill = "transparent", color = "black")+
  facet_wrap(~lyr, ncol = 3) +
  metR::scale_fill_divergent("m", na.value = "transparent")+
  guides(fill = guide_colourbar(theme = theme(
         legend.key.width  = unit(8, "lines"),
         legend.key.height = unit(0.5, "lines"),
         legend.title.position = "top",
         legend.title = element_text(hjust = 0.5)))) + 
  scale_x_continuous(limits = c(extent[1],extent[2]), expand = c(0,0))+
  scale_y_continuous(limits = c(extent[3],extent[4]), expand = c(0,0))+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")
        
```

```r
## Bootstrap heatwave tasmean reconstruction for a counterfactual/factual world
bs_sp_tx <- bs_spatanalogs(x = Tx,
                                     analogs = analogs$analogs_subperiods,
                                     n = 1000,
                                     event_fun = "mean", 
                                     detrend = F,
                                     replace = T,
                                     anom = T,
                                     ref_period = c(1950,2022))
                                     
counterfactual <- ((bs_sp_tx$`1951-1980`) %>% app("mean"))
factual <- ((bs_sp_tx$`1993-2022`) %>% app("mean"))
dif <- factual - counterfactual
all <- c(counterfactual,factual,dif)/10
names(all) <- c(names(bs_sp_tx),"dif")
extent <- ext(all)                

a_tx <- ggplot() +
  geom_spatraster(data = all[[-3]],interpolate = T) +
  geom_spatraster_contour_text(data = all[[-3]], breaks = seq(-1, 6, 1)) +
  geom_sf(data = borders, fill = "transparent", color = "black")+
  facet_wrap(~lyr, ncol = 3) +
  metR::scale_fill_divergent("ºC", na.value = "transparent")+
  guides(fill = guide_colourbar(theme = theme(
         legend.key.width  = unit(8, "lines"),
         legend.key.height = unit(0.5, "lines"),
         legend.title.position = "top",
         legend.title = element_text(hjust = 0.5)))) + 
  scale_x_continuous(limits = c(extent[1],extent[2]), expand = c(0,0))+
  scale_y_continuous(limits = c(extent[3],extent[4]), expand = c(0,0))+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")
        
b_tx <- ggplot() +
  geom_spatraster(data = all[[3]],interpolate = T) +
  geom_spatraster_contour_text(data = all[[3]], breaks = seq(-1, 6, 1)) +
  geom_sf(data = borders, fill = "transparent", color = "black")+
  facet_wrap(~lyr, ncol = 3) +
  scale_fill_steps2(midpoint = 0,low = "blue",mid = "white",high = "red",name = "ºC",
  na.value = "transparent", breaks = seq(-10,10,1), limits = c(0,4))+  
  guides(fill = guide_colourbar(theme = theme(
         legend.key.width  = unit(8, "lines"),
         legend.key.height = unit(0.5, "lines"),
         legend.title.position = "top",
         legend.title = element_text(hjust = 0.5)))) + 
  scale_x_continuous(limits = c(extent[1],extent[2]), expand = c(0,0))+
  scale_y_continuous(limits = c(extent[3],extent[4]), expand = c(0,0))+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")
       
# Z500 and T2m reconstructions composites 
(a_z500+b_z500 + plot_layout(widths = c(2, 1))) / (a_tx+b_tx + plot_layout(widths = c(2, 1)))
```


<img src="man/figures/bootstrap_maps_t2m_z500.png" alt="" width="700"/>

## Package citation {#package-citation}

...

## Contact {#contact}

Feel free to contact me: [marc.lemusicanovas\@eurac.edu](mailto:marc.lemusicanovas@eurac.edu){.email}
