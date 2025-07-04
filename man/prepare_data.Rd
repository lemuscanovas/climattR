% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_data.R
\name{prepare_data}
\alias{prepare_data}
\title{Prepare data from NetCDF with spatial and temporal processing}
\usage{
prepare_data(
  x,
  level = NULL,
  event_dates,
  time_window = 31,
  analog_months = NULL,
  detrend = FALSE,
  k = 2,
  scale = FALSE
)
}
\arguments{
\item{x}{A file path to a NetCDF file or a `SpatRaster` object.}

\item{level}{Optional; specify a particular pressure level to extract.}

\item{event_dates}{Dates of specific events to focus on within the data.}

\item{time_window}{Number of days before and after event_dates to include in the output.}

\item{analog_months}{Optional; specify months to be used in analysis if no specific time window provided.}

\item{detrend}{Logical; should the data be detrended? Defaults to FALSE.}

\item{k}{Degree of the polynomial for detrending if `detrend` is TRUE.}

\item{scale}{Logical; should the data be scaled? Defaults to FALSE.}
}
\value{
A list with two elements: `ts_wo_event`, the timeseries data excluding the event dates and window,
        and `event`, the timeseries data for the event dates.
}
\description{
This function processes NetCDF data with options for selecting specific pressure levels,
detrending, scaling, and extracting specific time windows around provided event dates.
It supports handling both daily and sub-daily resolutions.
}
\examples{
\dontrun{nc_file <- system.file("extdata", "example.nc", package = "yourPackageName")
event_dates <- as.Date(c("2022-06-01"))
result <- prepare_data(nc_file, event_dates = event_dates, time_window = 30)}
}
