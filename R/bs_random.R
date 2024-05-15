#' Bootstrap Randomization for Analog Subperiods
#'
#' This function performs bootstrap simulations based on subperiod analogs from time series data.
#' It offers options to detrend and normalize the data, convert hourly data to daily averages,
#' apply custom transformations, and compute statistical significance of the bootstrap results.
#'
#' @param x A data frame or tibble containing the time series data.
#' @param analogs A data frame of subperiods with associated analog information.
#' @param detrend Logical; if TRUE, the data is detrended using a polynomial of degree \code{k}.
#' @param k Integer; degree of the polynomial used for detrending (default is 2).
#' @param n Integer; the number of bootstrap samples to generate for each period.
#' @param replace Logical; if TRUE, sampling is done with replacement.
#' @param anom Logical; if TRUE, anomalies are calculated relative to a reference period.
#' @param ref_period Numeric vector; specifies the start and end years for the reference period if anomalies are to be calculated.
#' @param hour2day_fun String; specifies the function name to aggregate hourly data to daily (default is "mean").
#' @param event_fun Function; the function used to aggregate data during events.
#' @param conversion_fun Function; the function used to transform data before analysis.
#' @param cl Integer; the number of cores to use for parallel processing.
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item \code{observed}: The observed data filtered by the times of the analogs.
#'   \item \code{bootstrap_simulation}: A data frame of the bootstrap simulation results.
#'   \item \code{summary_bs}: A data frame summarizing the bootstrap results.
#' }
#'
#' @importFrom dplyr filter summarise mutate inner_join group_by ungroup select rename
#' @importFrom lubridate as_date ymd
#' @importFrom furrr future_map
#' @importFrom pbapply pblapply
#' @importFrom terra rast
#' 
#' @export

bs_random <- function(x, periods = c(1951,1980,1991,2020), n = 1000, event_fun = "mean", anom = F, 
                      ref_period = NULL, replace = T, detrend = F,k=2){
  
  event_FUN <- match.fun(FUN = event_fun)

  # Reading analogs dates and VOI nc ----------------------------------------
  dat <- x %>%
    rename("var" = 2)
  
  time_dat <- as_date(dat$time)
  
  
  if(is.Date(time_dat) == F){
    stop("time is not a date object. Please, check time values in your data")
  }
  
  # If hourly, converts to daily
  if(length(unique(dat$time)) != length(unique(time_dat))){
    dat <- dat %>%
      group_by(time) %>%
      summarise(var = mean(var, na.rm = T)) %>%
      ungroup()
  }    
  
  
  if(isTRUE(detrend)){
    detrend_pracma <- function(y, k) {
      fit <- polyfit(seq_along(y), y, k)
      y_pred <- polyval(fit, seq_along(y))
      y - y_pred
    }
    dat <- dat %>% na.omit() %>% mutate(var = detrend_pracma(var, k = k))
  }
  
  if(isTRUE(anom) & length(ref_period == 2)){
    ref_mean <- dat %>% 
      filter(year(time) >= ref_period[1], year(time) <= ref_period[2]) %>%
      summarise(ref_mean = mean(var, na.rm = T))
    
    dat$var <- dat$var - ref_mean$ref_mean
    
  }else if(isFALSE(anom)){
    
  } else if(isTRUE(anom) & is.null(ref_period)){
    stop("Please, a reference period is required. eg: ref_period = c(1981,2010)")
    
  }else{
    stop("Please, the reference period is as follows: ref_period = c(1981,2010)")
  }
  
  # Computing sd and mean for a counterfactual world -------------------------------
  
  # to select first period (factual)

  l = seq(1,length(periods),2)
  
  sim_bs_rnd_l <- list()
  for (ff in l) {
    
    dat_subset <- dat %>%
      filter(year(time) >= periods[ff],year(time) <= periods[ff+1]) %>%
      mutate(period = paste0(periods[ff],"-",periods[ff+1])) 
    
    bootstrap <- function(x){
      cases <- dat_subset %>%
        slice_sample(n = 1,replace = replace) %>%
        mutate(sim = paste0("sim",x))
    }
    
    message(paste0("bootstrapping the event ", n," times in the ", unique(dat_subset[,3]), ":"))
    
    sim_bs_rnd_l[[ff]] <- pbapply::pblapply(1:n, bootstrap, cl = 1) %>% 
      bind_rows()
    
  }
  sim_bs_rnd <- bind_rows(sim_bs_rnd_l)
  
  summary1_bootstrap_daily <- sim_bs_rnd %>%
    group_by(period) %>%
    summarise(var_mean = mean(var),
              var_median = median(var),
              var_sd = sd(var),
              .groups = "drop") %>%
    ungroup()
  
  return(list(bootstrap_simulation = sim_bs_rnd,
              summary_bs = summary1_bootstrap_daily))
  
}
