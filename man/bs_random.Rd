% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs_random.R
\name{bs_random}
\alias{bs_random}
\title{Bootstrap Randomization for Analog Subperiods}
\usage{
bs_random(
  x,
  periods = c(1951, 1980, 1991, 2020),
  n = 1000,
  event_fun = "mean",
  anom = F,
  ref_period = NULL,
  replace = T,
  detrend = F,
  k = 2
)
}
\arguments{
\item{x}{A data frame or tibble containing the time series data.}

\item{n}{Integer; the number of bootstrap samples to generate for each period.}

\item{event_fun}{Function; the function used to aggregate data during events.}

\item{anom}{Logical; if TRUE, anomalies are calculated relative to a reference period.}

\item{ref_period}{Numeric vector; specifies the start and end years for the reference period if anomalies are to be calculated.}

\item{replace}{Logical; if TRUE, sampling is done with replacement.}

\item{detrend}{Logical; if TRUE, the data is detrended using a polynomial of degree \code{k}.}

\item{k}{Integer; degree of the polynomial used for detrending (default is 2).}

\item{analogs}{A data frame of subperiods with associated analog information.}

\item{hour2day_fun}{String; specifies the function name to aggregate hourly data to daily (default is "mean").}

\item{conversion_fun}{Function; the function used to transform data before analysis.}

\item{cl}{Integer; the number of cores to use for parallel processing.}
}
\value{
A list containing three elements:
\itemize{
  \item \code{observed}: The observed data filtered by the times of the analogs.
  \item \code{bootstrap_simulation}: A data frame of the bootstrap simulation results.
  \item \code{summary_bs}: A data frame summarizing the bootstrap results.
}
}
\description{
This function performs bootstrap simulations based on subperiod analogs from time series data.
It offers options to detrend and normalize the data, convert hourly data to daily averages,
apply custom transformations, and compute statistical significance of the bootstrap results.
}
