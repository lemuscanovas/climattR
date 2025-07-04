#' Bootstrap Analogs for Time Series Data
#'
#' This function performs bootstrap simulations on time series data based on specified analog periods.
#' It includes options for detrending, calculating anomalies, and customizing bootstrap sampling.
#'
#' @param x A dataframe with at least 'time' and one variable column.
#' @param analogs A dataframe containing analogs information including 'time_obj' and 'period'.
#' @param n Integer; the number of bootstrap samples to generate for each period.
#' @param event_fun String; the name of the function to apply during event aggregation.
#' @param anom Logical; indicates whether anomalies relative to a reference period should be computed.
#' @param ref_period Numeric vector; the start and end years of the reference period for anomaly calculations.
#' @param replace Logical; indicates whether sampling should be with replacement.
#' @param detrend Logical; indicates whether the data should be detrended.
#' @param k Integer; the degree of polynomial used for detrending.
#'
#' @return A list containing three elements: observed data (aggregated original data), 
#'         bootstrap simulation results, and a summary of bootstrap results.
#'
#' @importFrom lubridate as_date
#' @importFrom terra rast
#' @importFrom dplyr filter summarise mutate inner_join group_by ungroup rename
#' @importFrom pbapply pblapply
#' @importFrom pracma polyfit polyval
#' @examples
#' \dontrun{
#' data <- data.frame(time = seq(as.Date('2020-01-01'), by = 'day', length.out = 365),
#'                    var = rnorm(365))
#' analogs <- data.frame(time_obj = seq(as.Date('2020-01-01'), by = 'month', length.out = 12),
#'                       period = rep('2020', 12))
#' results <- bs_analogs(data, analogs, n = 100, event_fun = "mean", anom = TRUE, ref_period = c(2019, 2021))
#' }
#' @export

bs_analogs <- function(x, 
                       analogs, 
                       n = 1000, 
                       event_fun = "mean", 
                       anom = F, 
                       ref_period = NULL,
                       replace = T, 
                       detrend = F, 
                       k = 2) {
  
  event_FUN <- match.fun(FUN = event_fun)

  # Reading analogs dates and VOI nc ----------------------------------------
  dat <- x %>%
    rename("var" = 2)

  time_dat <- as_date(dat$time)
  
  
  if(is.Date(time_dat) == F){
    stop("time is not a date object. Please, check time values in your data")
  }
  
  # If hourly, converts to daily
  message("Caution! if hourly data is provided then it is converted into daily mean.")
  
  if(length(unique(dat$time)) != length(unique(time_dat))){
      dat <- dat %>%
        group_by(time) %>%
        summarise(var = mean(var)) %>%
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
  
    }else if(isFALSE(anom)){
      
    } else if(isTRUE(anom) & is.null(ref_period)){
    stop("Please, a reference period is required. eg: ref_period = c(1981,2010)")
    
  }else{
    stop("Please, the reference period is as follows: ref_period = c(1981,2010)")
      }

  # Computing sd and mean for a counterfactual world -------------------------------
  
  # to select first period (factual)
  yr_split <- analogs$period %>% 
    unique() 
  
  # selecting n analogs per period 
  sim_bs_l <- list()
  for(ii in seq_along(yr_split)){
    
  dat_subset <- inner_join(analogs, dat, by = "time") %>%
    filter(str_detect(period, (yr_split %>% str_sub(1,4))[ii])) %>% 
    relocate(period, .after = time)
  
  bootstrap <- function(x){
    cases <- dat_subset %>% group_by(time_obj) %>%
      slice_sample(n = 1,replace = replace) %>%
      mutate(sim = paste0("sim",x))
  }
  
  message(paste0("bootstrapping the event ", n," times in the ", yr_split[ii], ":"))
  
  sim_bs_l[[ii]] <- pbapply::pblapply(1:n, bootstrap, cl = 1) %>% 
                bind_rows()

  }
  
  sim_bs <- bind_rows(sim_bs_l)
 
  
  # COMPUTING STATISTICAL SIGN OF BOOTSTRAP  --------------------------------
  
  .range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  
  
  join_bootstraps <- sim_bs %>% 
    dplyr::select(sim,time_obj, dist, var,period) %>%
    # group_by(time_obj) %>%
    # mutate(dist = .range01(dist,na.rm = T)) %>%
    group_by(time_obj,sim,period) %>%
    summarise(var = ifelse(exists("ref_mean"), var - pull(ref_mean,1),var),
           dist = mean(dist),
           var = var,.groups = "drop") %>%
    ungroup()
  
  # mean and sd of RMSD or Euclid dist, 
  # mean, median and sd of var
  summary1_bootstrap_daily <- join_bootstraps %>%
    group_by(time_obj,period) %>%
    summarise(dist_mean = mean(dist),
              dist_sd = sd(dist),
              var_mean = mean(var),
              var_median = median(var),
              var_sd = sd(var),
              .groups = "drop") %>%
    ungroup()
  
  # mitjana anomalia de cada dia analeg
  # ttest of metric distance and var
  # summary2_bootstrap_daily <- join_bootstraps %>%
  #   pivot_wider(
  #     names_from = period,
  #     values_from = c(dist, var),values_fn = list) %>%
  #   unnest(cols = 2:5) %>%
  #   group_by(time_obj) %>%
  #   summarise(ttest_dist = t.test(.[[2]],.[[3]])$p.value,
  #             ttest_mean = t.test(.[[4]],.[[5]])$p.value)
  # 
  # summary2_bootstrap_event <- join_bootstraps %>%
  #   pivot_wider(
  #     names_from = period,
  #     values_from = c(dist, var),values_fn = list) %>%
  #   unnest(cols = 2:5) %>%
  #   summarise(ttest_dist = t.test(.[[2]],.[[3]])$p.value,
  #             ttest_mean = t.test(.[[4]],.[[5]])$p.value)
  

  dat_days_event <- dat %>% 
    filter(time %in% unique(analogs$time_obj)) %>%
    group_by(time) %>%
    summarise(var = ifelse(exists("ref_mean"), var - pull(ref_mean,1),var),.groups = "drop") %>%
    ungroup() %>%
    rename(time_obj = time) 
  

  return(list(observed = dat_days_event,
              bootstrap_simulation = join_bootstraps,
              summary_bs = summary1_bootstrap_daily))
         
}
