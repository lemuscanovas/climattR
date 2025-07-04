% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analogs_searcher.R
\name{analogs_searcher}
\alias{analogs_searcher}
\title{Search for analogs based on historical data and specific events}
\usage{
analogs_searcher(
  ts_wo_event,
  event,
  n = 20,
  periods = c(1951, 1980, 1991, 2020),
  metric = "rmsd"
)
}
\arguments{
\item{ts_wo_event}{List of `SpatRaster` objects containing timeseries data excluding the event dates.}

\item{event}{List of `SpatRaster` objects containing timeseries data for the event dates.}

\item{n}{The number of top analogs to return.}

\item{periods}{A vector of years to define the periods for analysis, should be in pairs.}

\item{metric}{A character string specifying the metric for comparison: 'rmsd', 'euclidean', or 'pearson'.}
}
\value{
A list containing the analogs for the full period and, if specified, for subperiods.
}
\description{
This function compares historical meteorological data with event data to find analogs.
It supports multiple metrics for comparison including RMSD, Euclidean, and Pearson correlation.
Analog dates are identified based on the smallest distance metrics calculated over specified periods.
}
\examples{
\dontrun{ts_wo_event_example <- list(rast(system.file("extdata", "example1.nc", package = "yourPackageName")))
event_example <- list(rast(system.file("extdata", "example2.nc", package = "yourPackageName")))
result <- analogs_searcher(ts_wo_event_example, event_example, n = 5, periods = c(1951,1980,1991,2020), metric = "rmsd")
}
}
